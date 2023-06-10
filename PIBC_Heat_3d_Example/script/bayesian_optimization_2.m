function sol = bayesian_optimization_2(data,N_initial,N_evaluate,true_theta)
d = data.d*(data.d+1)/2;
N_initial = N_initial * d;
N_evaluate = N_evaluate*d;
if nargin == 3
    true_theta = parametrization(data.D_true);
end
D_theta = data.design_para(1:N_initial,:);
obj = @(x) mse_pde(x, data);
results = my_EI(obj,D_theta,N_evaluate-N_initial,data);
%[x2,fval2,exitflag2] = patternsearch(@f_opt_emu,[x01],[],[],[],[],range(1,:),range(2,:),[],options_patternsearch);
sol.theta_opt = results.x_opt;
sol.theta_true = true_theta;

sol.theta_opt_trace = results.X_trace;
sol.y_trace = results.Y_trace;
sol.RMSE = mean((true_theta-sol.theta_opt_trace).^2,2).^(1/2)';
sol.sigmasq = data.sigmasq;
end


function result = my_EI(obj,D_ini,N_evaluate,data)
d = data.d*(data.d+1)/2;
N_initial = size(D_ini,1);
for i = 1:N_initial
simu(i) = obj(D_ini(i,:));
end
X_trace = zeros(N_evaluate+N_initial, d);
Y_trace = zeros(1, N_evaluate+N_initial);
x_current = D_ini (simu == min(simu),:);
y_best = min(simu);
GP_model = fit_GP(D_ini,simu);
Y_trace(1:N_initial) = ones(N_initial,1) * y_best;
X_trace(1:N_initial,:) = ones(N_initial,1) * x_current;
for i =1 : N_evaluate
    onj_EI_opt = @(x) obj_EI(x, y_best, GP_model);
    %x_next = fmincon(onj_EI_opt, x_current, [],[],[],[],data.lb,data.ub);
    a = haltonset(d);
    a = net(a,5*d);
    start_set = a.*(data.ub-data.lb-0.02)+data.lb+0.01;
    problem = createOptimProblem('fmincon','x0',start_set(1,:), 'lb',data.lb,'ub',data.ub,'objective',onj_EI_opt);
    tpoints=CustomStartPointSet(start_set);
    ms=MultiStart('StartPointsToRun','all','Display','off');
    [x01,fva10,exitflag0]=run(ms,problem,tpoints);
    options = optimoptions('particleswarm','Display','off','MaxIterations',1000,'HybridFcn',@patternsearch);
    [x_next,fval1,exitflag0]= particleswarm(onj_EI_opt,d,data.lb,data.ub,options);
    if fva10<fval1
        x_next = x01;
    end
    y_next = mse_pde(x_next, data);
    %[x_next,y_next]
    simu = [simu,y_next];
    D_ini = [D_ini; x_next];
    GP_model = fit_GP(D_ini,simu);
    x_current = D_ini (simu == min(simu),:);
    y_best = min(simu);
    X_trace(i+N_initial,:) = x_current;
    %[a, b] = predict(GP_model,x_current)
    Y_trace(i+N_initial) = y_best;
end
result.x_opt = x_current;
result.y_opt = y_best;
result.X_trace = X_trace;
result.Y_trace = Y_trace;
end


function y = obj_EI(x_pred,y_best,GP_model)
[pred_mean, pred_sd] = predict_GP(GP_model,x_pred);
stand_err = (y_best - pred_mean) / pred_sd;
y = -(y_best - pred_mean) * normcdf(stand_err) - pred_sd * normpdf(stand_err); 
end


function GP_model = fit_GP(X,Y)
d = size(X,2);
Dist = cell(1,d);
for i =1:d 
    Dist{i} = abs(X(:,i)-X(:,i)');
end
obj = @(theta)likelihood(theta,Dist,Y);
scale=[range(X,1)];
options=optimoptions(@fminunc,'MaxIter',10^3,'TolX',0.01,'TolFun',0.01,'MaxFunEvals',10^5,'Display','off','Algorithm','quasi-newton','ObjectiveLimit',-10^250);
[paropt1, fval1, exitflag1]=fminunc(obj,[log(0.909700041540068.*scale)],options);  
[paropt2, fval2, exitflag2]=fminunc(obj,[log(2.54815755509488.*scale)],options);   
if fval1<fval2
    theta_opt = paropt1;
else 
    theta_opt = paropt2;
end


Corr = corr_matrix(theta_opt,Dist);
L = chol(Corr+1e-6*eye(size(Corr,1)));
Y_t = (Y-mean(Y))*inv(L);
GP_model.R = Corr;
GP_model.sigma = mean(Y_t.^2);
GP_model.Y = Y;
GP_model.X = X;
GP_model.theta = theta_opt;

end


function y = likelihood(theta,Dist,Y)
Corr = corr_matrix(theta,Dist);
[L,flag] = chol(Corr+1e-6*eye(size(Corr,1)));
Y_t = (Y-mean(Y))/inv(L);
y = sum(log(diag(L))) + sum(Y_t.^2)/2;
end


function Corr = corr_matrix(theta,Dist)
d = length(Dist);
theta = exp(theta);
Corr = 1;
for i = 1:d
    rho1=sqrt(6)*Dist{i}./theta(i);
    Corr=Corr.*((exp(-rho1)).*(rho1+1));
end
end

function [y_new,sd_new] = predict_GP(GP_model,x_new)
d = size(GP_model.X,2);
Dist = cell(1,d);
for i =1:d 
    Dist{i} = abs(x_new(:,i)-GP_model.X(:,i)');
end
dist = pdist2(x_new,GP_model.X);
r = corr_matrix(GP_model.theta,Dist) + 1e-6 *(dist ==0);
R = GP_model.R + 1e-6 * eye(size(GP_model.R,1));
Y = GP_model.Y;
y_new = mean(Y)+ r* (R\(Y'-mean(Y)));
sd_new = GP_model.sigma * (1-r*(R\r'));
sd_new = (max(0,sd_new))^(1/2);
end







