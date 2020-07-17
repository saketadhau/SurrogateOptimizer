
import casadi.*

clear
clc
global nx nu nd
global lbx ubx dx0 lbu ubu u0
par.tf = 1/60;

% parameters
par.Cs = 0.0549 ; % Heat capacities
par.Ci = 0.0928;
par.Ce = 3.32;
par.Ch = 0.889;

par.Ris = 1.89;
par.Rih = 0.146;
par.Rie = 0.897;
par.Rea = 4.38;
par.Ria = 2.5;

par.Aw = 5.75;
par.Ae = 3.87;

[sys,par] = BuildingModel(par);

nx = numel(sys.x);
nu = numel(sys.u);
nd = numel(sys.d);

lbx = -60.*ones(nx,1);
ubx = 100*ones(nx,1);
dx0 = 20*ones(nx,1);
lbu = 0*ones(nu,1);
ubu = 100*ones(nu,1);
u0  = 0;

u_in = 8;
d_val = [0.0,1,22]';
xf = solveODE(sys,d_val,u_in);


NN = 1;

if ~NN
    par.N = 3/par.tf;
    [solver,par] = buildNLP(sys.f,par);
end

par.nIter = 12*60;

for sim_k = 1:par.nIter
    
    if sim_k>3*60 && sim_k<=12*60
        d_val(3) = 22;
    else
        d_val(3) =18;
    end
    
    if sim_k>9*60
        d_val(2)= 4;
    end
    
    if sim_k>6*60
        d_val(1)= 0.1;
    end
    NMPC.d(:,sim_k) = d_val;
    
    
    
    if NN
        tic;
        NMPC.u(:,sim_k) = ApproxMPC(vertcat(xf,u_in,d_val));
        par.elapsednlp = toc;
        
        u_in = NMPC.u(:,sim_k);
        NMPC.sol_t(sim_k) = par.elapsednlp;
    else
        tic;
        sol = solver('x0',par.w0,'p',vertcat(xf,u_in,d_val),'lbx',par.lbw,'ubx',par.ubw,'lbg',par.lbg,'ubg',par.ubg);
        par.elapsednlp = toc;
        
        flag = solver.stats();
        exitflag =  flag.return_status;
        disp(['Data point nr. ' num2str(sim_k)  ' - ' exitflag '. CPU time: ' num2str(par.elapsednlp) 's'])
        
        w_opt = full(sol.x);
        
        NMPC.sol_t(sim_k) = par.elapsednlp;
        
        n_w_i = nx + par.N*(4*nx+nu);
        
        u1_opt = [w_opt(nx+1:4*nx+nu:n_w_i);NaN];
        x1_opt = w_opt([1,nu+4*nx+1:4*nx+nu:n_w_i]);
        x2_opt = w_opt([2,nu+4*nx+2:4*nx+nu:n_w_i]);
        x3_opt = w_opt([3,nu+4*nx+3:4*nx+nu:n_w_i]);
        x4_opt = w_opt([4,nu+4*nx+4:4*nx+nu:n_w_i]);
        
        NMPC.u(:,sim_k) = [w_opt(nx+1)];
        u_in = NMPC.u(:,sim_k);
        
    end
    %------------------------ Plant simulation----------------------
    
    Fk = sys.F('x0',xf,'p',vertcat(u_in,d_val));
    xf =  full(Fk.xf) ;
    
    NMPC.J(sim_k) = full(Fk.qf);
    NMPC.x(:,sim_k) = xf;
    
    
end


%% plotting

figure(56)
clf
subplot(211)
hold all
plot(NMPC.d(3,:),'k','linewidth',2)
plot(NMPC.x(2,:),'linewidth',2)
ylabel('$T_i$','Interpreter','latex')
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
box on

subplot(212)
hold all
stairs(NMPC.u,'linewidth',2)
ylabel('$u$','Interpreter','latex')
axs = gca;
axs.FontSize = 14;
axs.TickLabelInterpreter = 'latex';
box on


