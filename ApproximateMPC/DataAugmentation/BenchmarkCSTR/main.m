import casadi.*

clear
clc

global nx nu nd
global lbx ubx dx0 lbu ubu u0
par.tf = 3;

[sys,par] = BenchmarkCSTR(par);

nx = numel(sys.x);
nu = numel(sys.u);
nd = numel(sys.d);


lbx = 0.*ones(nx,1);
ubx = 1*ones(nx,1);
dx0 = [0.98;0.39];
lbu = 0*ones(nu,1);
ubu = 2*ones(nu,1);
u0  = 1;

u_in = 0.1;
d_val = [0.2632;0.6519];
xf = [0.2;0.4];

NN = 1;

if ~NN
    par.N = 140;
    [solver,par] = buildNLP(sys.f,par);
end

par.nIter = 140;

for sim_k = 1:par.nIter
    
    
    if NN
        % Approximate MPC
        tic;
        NMPC.u(:,sim_k) = max(0,min(2,ApproxMPC(vertcat(xf,d_val))));
        par.elapsednlp = toc;
        
        u_in = NMPC.u(:,sim_k);
        NMPC.sol_t(sim_k) = par.elapsednlp;
    else
        % Traditional MPC
        tic;
        sol = solver('x0',par.w0,'p',vertcat(xf,d_val),'lbx',par.lbw,'ubx',par.ubw,'lbg',par.lbg,'ubg',par.ubg);
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
        
        NMPC.u(:,sim_k) = [w_opt(nx+1)];
        u_in = NMPC.u(:,sim_k);
        
    end
    
    %------------------------ Plant simulation----------------------
    
    Fk = sys.F('x0',xf,'p',vertcat(u_in,d_val));
    xf =  full(Fk.xf) ;
    
    NMPC.J(sim_k) = full(Fk.qf);
    NMPC.x(:,sim_k) = xf;
    
    
end

save('nmpc','NMPC')

figure(56)
hold all
stairs(NMPC.u,'linewidth',2)
ylabel('$u$','interpreter','latex')
axs = gca;
axs.FontSize = 14
axs.TickLabelInterpreter = 'latex';
box on




