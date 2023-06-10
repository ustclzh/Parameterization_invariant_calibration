function [D,x2,hessian]=map_optimization(data,lik)
    if nargin ==1
        global range data_opt
        range=[0.01,-1,0.01,-1,-1,0.01;1,1,1,1,1,1]*10^(1/2);
        d=data.d;
        range=range(:,1:d*(d+1)/2);
        data_opt=data;
        options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000);
        options = optimoptions('particleswarm','Display','off','MaxIterations',3000,'HybridFcn',@patternsearch);
        [x01,fval0,exitflag0]= particleswarm(@f_opt_emu,d*(d+1)/2,range(1,:),range(2,:),options);
        [x2,fval2,exitflag2] = patternsearch(@f_opt_emu,[x01],[],[],[],[],range(1,:),range(2,:),[],options_patternsearch);
        options_neal=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter', 5000);
        [x3,fval3]=fminsearch(@f_opt_emu,x2,options_neal);
        options_bnd=optimset('TolX',1e-6,'MaxIter', 5000);
        [x4,fval4]=fminbnd(@f_opt_emu,-4,4,options_bnd);

        [fval0,fval2,fval3,fval4]
    %     n=size(data.design_matrix,3);
    %     start_set=[data.design_para((n-10):n,:);x2]
    %     problem = createOptimProblem('fminunc','x0',ones(1,d*(d+1)/2),'objective',@f_opt_emu);
    %     tpoints=CustomStartPointSet(start_set(:,1:d*(d+1)/2));
    %     ms=MultiStart('StartPointsToRun','all');
    %     [x02,fval0,exitflag0]=run(ms,problem,tpoints);
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
    else
        global range data_opt
        range=[0.1,-1,0.1,-1,-1,0.1;1,1,1,1,1,1]*10^(1/2);
        d=data.d;
        range=range(:,1:d*(d+1)/2);
        data_opt=data;
        options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000);
        options = optimoptions('particleswarm','Display','off','MaxIterations',3000,'HybridFcn',@patternsearch);
        [x01,fval0,exitflag0]= particleswarm(@f_opt_emu_lik,d*(d+1)/2,range(1,:),range(2,:),options);
        [x2,fval2,exitflag2] = patternsearch(@f_opt_emu_lik,[x01],[],[],[],[],range(1,:),range(2,:),[],options_patternsearch);
    %     n=size(data.design_matrix,3);
    %     start_set=[data.design_para((n-10):n,:);x2]
    %     problem = createOptimProblem('fminunc','x0',ones(1,d*(d+1)/2),'objective',@f_opt_emu);
    %     tpoints=CustomStartPointSet(start_set(:,1:d*(d+1)/2));
    %     ms=MultiStart('StartPointsToRun','all');
    %     [x02,fval0,exitflag0]=run(ms,problem,tpoints);
        D=invparametrization(x2,data_opt.parameterization);
        if max(eig(D))>=10 || min(eig(D))<=0.1
            D0=5*eye(d);
            x2p=parametrization(D0,data_opt.parameterization);
            [x0,fval,exitflag,output,grad,hessian] = fminunc(@f_opt_emu1_lik,x2p);
        else
            [x0,fval,exitflag,output,grad,hessian] = fminunc(@f_opt_emu1_lik,x2);
        end
        
        if fval<fval2
            D=invparametrization(x0,data_opt.parameterization);
            x2=x0;
        end

    end
end


function y=f_opt_emu(L)
    global data_opt
    D=invparametrization(L,data_opt.parameterization);
    y=-log_posterior_emulator(data_opt,D);
    e=eig(D);
    y=y + 10^6*(e(1)<0.1) + 10^6*(e(end)>10);
end

function y=f_opt_emu1(L)
    global data_opt
    D=invparametrization(L,data_opt.parameterization);
    if max(eig(D))>=10 || min(eig(D))<=0.1
        y=inf;
        return
    end
    
    y=-log_posterior_emulator(data_opt,D);
end


function y=f_opt_emu_lik(L)
    global data_opt
    D=invparametrization(L,data_opt.parameterization);
    y=-log_likelihood_emulator(data_opt,D);
    e=eig(D);
    y=y + det(D)^(-(size(D,1)+1)/2)+10^6*(e(1)<0.1) + 10^6*(e(end)>10);
end

function y=f_opt_emu1_lik(L)
    global data_opt
    D=invparametrization(L,data_opt.parameterization);
    if max(eig(D))>=10 || min(eig(D))<=0.1
        y=inf;
        return
    end
    y=-log_likelihood_emulator(data_opt,D);
    y=y+det(D)^(-(size(D,1)+1)/2);
end