function [D,x2,hessian]=mle_optimization(data)
    global range data_opt
    range=[0.1,-1,0.1,-1,-1,0.1;1,1,1,1,1,1]*10^(1/2);
    d=data.d;
    range=range(:,1:d*(d+1)/2);
    data_opt=data;
    options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000);
    options = optimoptions('particleswarm','Display','off','MaxIterations',3000,'HybridFcn',@patternsearch);
    [x01,fval0,exitflag0]= particleswarm(@f_opt_emu,d*(d+1)/2,range(1,:),range(2,:),options);
    [x2,fval2,exitflag2] = patternsearch(@f_opt_emu,[x01],[],[],[],[],range(1,:),range(2,:),[],options_patternsearch);
    D=invparametrization(x2,data_opt.parameterization);
    if max(eig(D))>10 || min(eig(D))<0.1
        D0=5*eye(d);
        x2p=parametrization(D0,data_opt.parameterization);
        [x0,fval,exitflag,output,grad,hessian] = fminunc(@f_opt_emu1,x2p);
    else
        [x0,fval,exitflag,output,grad,hessian] = fminunc(@f_opt_emu1,x2);
    end
    if fval<fval2
        D=invparametrization(x0,data_opt.parameterization);
        x2=x0;
    end
end


function y=f_opt_emu(L)
    global data_opt
    D=invparametrization(L,data_opt.parameterization);
    y=-log_likelihood_emulator(data_opt,D);
end


function y=f_opt_emu1(L)
    global data_opt range
    if sum(L<range(1,:))>0
        L=range(1,1:length(L));
    end
    D=invparametrization(L,data_opt.parameterization);
    y=-log_likelihood_emulator(data_opt,D);
end
