function data=config(data,d,large_initial)


data.parameterization=1; % parametrizatipn of tensor, 1 for cholesky, 2 for matrix logrithm
data.distance=1; % distance matrix used to construct kernel of GP, 1 for 
data.iso_kernel=1; %1 isotopic kernel, 0 product kernel; 
data.cor_func=1; % correlation funciton, 1 for 3/2 matern
data.design_algo=1; % design algorithm. 1 for support point, 2 for finite region
if nargin==2
    data.q_0=5*d*(d+1)/2; % size of initial design
else
    data.q_0=10*d*(d+1)/2; % size of initial design
end
data.q_1=5*d*(d+1)/2; % size of follow up design
data.eig_L=0.01;
data.eig_U=10;
data.d=d;
data.N_MCMC=20000;
end