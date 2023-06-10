function [D,x3,hessian]=map_optimization(data,strategy,V,v)
switch strategy 
    case 1 % emulator based 
        f_obj = @(L) f_posterior_emu(data, L, V, v);
    case 2 % simulator based 
        f_obj = @(L) f_posterior_simu(data, L, V, v);
    case 3 % likelihood emulator
        f_obj = @(L) f_likelihood_emu(data, L);
    case 4 % likelihood simulator
        f_obj = @(L) f_likelihood_simu(data, L);
end

d=data.d;
range=design_region(data);
        range=[0.1,-1,0.1,-1,-1,0.1;1,1,1,1,1,1]*10^(1/2);
        d=data.d;
        range=range(:,1:d*(d+1)/2);
options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000);
options = optimoptions('particleswarm','Display','off','MaxIterations',3000,'HybridFcn',@patternsearch);
[x01,fval0,exitflag0]= particleswarm(f_obj,d*(d+1)/2,range(1,:),range(2,:),options);
[x2,fval2,exitflag2] = patternsearch(f_obj,[x01],[],[],[],[],range(1,:),range(2,:),[],options_patternsearch);
options_neal=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter', 5000, 'Display', 'off');
[x3,fval3]=fminsearch(f_obj,x2,options_neal);

L_true = parameterization(data.D_true);
[x4,fval4]=fminsearch(f_obj,L_true,options_neal);
if fval4 <=fval3
x3 = x4;
fval3 = fval4;
end
%options_bnd=optimset('TolX',1e-6,'MaxIter', 5000);
%[x4,fval4]=fminbnd(f_obj,-4,4,options_bnd);
%[exitflag0,exitflag2]

n=size(data.design_matrix,3);
if n<13
    start_set=[data.design_para;x3;L_true];
else
    start_set=[data.design_para((n-12):n,:);x3;L_true];
end
problem = createOptimProblem('fminunc','x0',ones(1,d*(d+1)/2),'objective', f_obj);
tpoints=CustomStartPointSet(start_set(:,1:d*(d+1)/2));
ms=MultiStart('StartPointsToRun','all','Display','off');
[x02,fval4,exitflag0]=run(ms,problem,tpoints);
%[fval0,fval2,fval3,fval4]
if fval4 <fval3
    x3 = x02;
    fval3 = fval4;
end

D=invparametrization(x3,data.parameterization);
option_fminunc = optimoptions('fminunc', 'Display', 'off');
if max(eig(D))>100 || min(eig(D))<0.01
    D0=eye(d);
    x2p=parametrization(D0,data.parameterization);
    [x0,fval,exitflag,output,grad,hessian] = fminunc(f_obj,x2p,option_fminunc);
else
    [x0,fval,exitflag,output,grad,hessian] = fminunc(f_obj,x3,option_fminunc);
end

if fval<fval3
    D=invparametrization(x0,data.parameterization);
    x3=x0;
end
end

function y=f_posterior_emu(data, L, V, v)
D=invparametrization(L, data.parameterization);
if ispd(data, D) == 0
    y = 10^256;
return
end
y=-log_posterior_emulator(data,D, V, v);
end


function y=f_likelihood_emu(data, L)
D=invparametrization(L,data.parameterization);
if ispd(data, D) == 0
    y = 10^256;
return
end
y=-log_likelihood_emulator(data,D);
end

function y=f_posterior_simu(data, L, V, v)
D=invparametrization(L, data.parameterization);
if ispd(data, D) == 0
    y = 10^256;
    return;
end
if nargin == 2
    y=-log_posterior_simulator(data,D);
else
    y=-log_posterior_simulator(data,D,V, v);
end
end

function y=f_likelihood_simu(data, L)
D=invparametrization(L, data.parameterization);
if ispd(data, D) == 0
    y = 10^256;
    return;
end
y=-log_likelihood_simulator(data,D);
end

function y = ispd(data, D)
e = eig(D);
if min(e)>=data.eig_L && max(e)<=data.eig_U
    if issymmetric(D)
        y = 1;
        return
    end
end
y = 0;
end