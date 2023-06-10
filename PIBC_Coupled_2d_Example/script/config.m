function data=config(eg,distance_gp,initial_size,sequential_size)

% setup simulator information
data.eg = eg;

% setup for parameter 
d = 2 + (eg == 3) + (eg == 4)+ (eg == 6)+ (eg == 7) + (eg == 9);
data.d=d;

data.parameterization=1; % parametrizatipn of tensor, 1 for cholesky, 2 for matrix logrithm
data.distance=1; % distance matrix used to construct kernel of GP, 1 for 
data.iso_kernel=1; %1 isotopic kernel, 0 product kernel; 
data.cor_func=1; % correlation funciton, 1 for 3/2 matern
data.design_algo=2; % design algorithm. 1 for support point, 2 for finite region
data.q_0=initial_size*d*(d+1)/2; % size of initial design
data.q_1=sequential_size*d*(d+1)/2; % size of follow up design
data.distance_gp = distance_gp; % use how many decomposition

% setup for MCMC
data.N_MCMC=20000;

% boundary of the 
data.eig_L=0.1;
data.eig_U=20;
data.lb = [data.eig_L,-data.eig_U,data.eig_L,-data.eig_U,-data.eig_U,data.eig_L];
data.ub = data.eig_U*ones(1,6);
data.theta_bo=[
        optimizableVariable('theta1', [0.01,data.eig_U], 'type', 'real')
        optimizableVariable('theta2', [-data.eig_U,data.eig_U], 'type', 'real')
        optimizableVariable('theta3', [0.01,data.eig_U], 'type', 'real')
        optimizableVariable('theta4', [-data.eig_U,data.eig_U], 'type', 'real')
        optimizableVariable('theta5', [-data.eig_U,data.eig_U], 'type', 'real')
        optimizableVariable('theta6', [-data.eig_U,data.eig_U], 'type', 'real')];
if d ==2
    data.lb = data.lb(1:3);
    data.ub = data.ub(1:3);
    data.theta_bo=[optimizableVariable('theta1', [0.01,data.eig_U], 'type', 'real')
        optimizableVariable('theta2', [-data.eig_U,data.eig_U], 'type', 'real')
        optimizableVariable('theta3', [0.01,data.eig_U], 'type', 'real')];
end


end