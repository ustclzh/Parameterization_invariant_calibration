function sol = bayesian_optimization(data,N_initial, N_evaluate,true_theta)
d = data.d*(data.d+1)/2;
N_initial = N_initial * d;
N_evaluate = N_evaluate * d;
if nargin == 3
    true_theta = parametrization(data.D_true);
end
theta_bo =data.theta_bo;
obj = @(x) mse_pde_bo(x,data);
results = bayesopt(obj,theta_bo,'MaxObjectiveEvaluations',N_evaluate, 'NumSeedPoints',N_initial,'Verbose',1,'PlotFcn',{});
X_trace= table2array(results.XTrace);
X_opt_trace = X_trace(results.IndexOfMinimumTrace,:);
theta_opt =table2array(results.XAtMinObjective);


sol.theta_opt = theta_opt;
sol.theta_true = true_theta;
sol.y_trace = results.ObjectiveMinimumTrace';
sol.theta_opt_trace = X_opt_trace;
sol.RMSE = mean((true_theta-sol.theta_opt_trace).^2,2).^(1/2)';
sol.sigmasq = data.sigmasq;
end


function loss = mse_pde_bo(theta, data)
eg = data.eg;
if eg ==1
    theta = [theta.theta1,theta.theta2,theta.theta3];
elseif eg==2
    theta = [theta.theta1,theta.theta2,theta.theta3];
elseif eg ==3
    theta = [theta.theta1,theta.theta2,theta.theta3,theta.theta4,theta.theta5,theta.theta6];
elseif eg ==4 
    theta = [theta.theta1,theta.theta2,theta.theta3,theta.theta4,theta.theta5,theta.theta6];
elseif eg ==5 
    theta = [theta.theta1,theta.theta2,theta.theta3];
elseif eg ==6 
    theta = [theta.theta1,theta.theta2,theta.theta3,theta.theta4,theta.theta5,theta.theta6];
elseif eg ==7
    theta = [theta.theta1,theta.theta2,theta.theta3,theta.theta4,theta.theta5,theta.theta6];
end
loss = mse_pde(theta, data);
end




