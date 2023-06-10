function RESULT=run_simulator_based_infer(data,mcmc)
global range data_opt
range=[0.01,-1,0.01,-1,-1,0.01;1,1,1,1,1,1]*10^(1/2);
d=data.d;
range=range(:,1:d*(d+1)/2);
data_opt=data;
options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000);
options = optimoptions('particleswarm','Display','off','MaxIterations',3000,'HybridFcn',@patternsearch);
[x01,fval0,exitflag0]= particleswarm(@f_opt,d*(d+1)/2,range(1,:),range(2,:),options);
[x2,fval2,exitflag2] = patternsearch(@f_opt,[x01],[],[],[],[],range(1,:),range(2,:),[],options_patternsearch);
[x,fval,exitflag,output,grad,hessian] = fminunc(@f_opt,x2);
disp('MAP of Simulator based Posterior')
D=invparametrization(x,data_opt.parameterization)
% MCMC based on simulator
if nargin == 1
    para_MCMC.N=10000;
    start_point=parametrization(D,data_opt.parameterization);
    result_simulator=MCMC_tensor(data_opt,para_MCMC,start_point);

    RESULT.result_simulator=result_simulator;
end
RESULT.D=D;
RESULT.L = x;
RESULT.hessian=hessian;
end

function y=f_opt(L)
global data_opt
D=invparametrization(L,data_opt.parameterization);
y=-log_posterior_simulator(data_opt,D);
end