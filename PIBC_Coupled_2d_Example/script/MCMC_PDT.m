function result=MCMC_PDT(data,para_MCMC,start_point,Sigma, prior_strategy)
    global MCMC_scale
    MCMC_scale=2.38/(4*data.d*(data.d+1))^(1/2);
    N=para_MCMC.N;
    v_hyper = data.d+1;
    V0 = eye(data.d);
    %start_point=data.design_para(end,:);
    V=eye(data.d);
    v0=data.d+1;
    M=length(start_point);
    chain=zeros(N,M);
    chain_V=cell(1,N);
    chain_V_L=zeros(N,M);
    chain_v0=zeros(N,1);
    L_current=start_point;
    D_current=invparametrization(L_current,data.parameterization);
    %logpost_current=log_posterior_emulator(data,D_current,V,v0);
    if nargin==3 % specify stating point
        [y,var]=Emulator_Predict(data,D_current);
        Grad=gradient_tensor(data,L_current);
        n=size(Grad,2);
        k=size(L_current,2);
        data.Sigma=invandlogdet(Grad*((data.sigmasq*eye(n)+var)\Grad')+eye(k));
        prior_strategy = 1;
    else % specify the starting point for 
        data.Sigma=Sigma;
    end
    if prior_strategy == 0
        v0 = 0;
    elseif prior_strategy == 2
        v0 = data.d+1;
    end
    accep=0;
    fprintf('\n')
    time_used=0;
    D=invparametrization(L_current,data.parameterization);
    for i=1:N
        tic
        logpost_current=log_posterior_emulator(data,D_current,V,v0);
        [L_new,Sigma,Sigma_new]=proposal(data,L_current,V,size(V,1));
        %L_new=proposal(data,L_current);
        D_new=invparametrization(L_new,data.parameterization);
        logpost_new=log_posterior_emulator(data,D_new,V,v0);
        alpha=exp(logpost_new-logpost_current)*mvnpdf(L_current,L_new,Sigma_new)/mvnpdf(L_new,L_current,Sigma);
        u=rand(1);
        if chol_bound(D_new)==1
            chain(i,:)=L_current;
        elseif alpha>u % accept
            chain(i,:)=L_new;
            L_current=L_new;
            D=D_new;
            accep=accep+1;
            logpost_current=logpost_new;
        else
            chain(i,:)=L_current;
        end
        time_Z=toc;
        time_used=time_used+time_Z;
        if rem(i,100)==0
            fprintf(repmat('\b',1,200))
            fprintf('MCMC %5d 00/ %5d 00 Steps, Time: %3.3f, Accept Rate: %2.2f',floor(i/100), floor(N/100) ,time_used,100*accep/i)
        end
        if prior_strategy == 1
            V=iwishrnd(D_current+V0,v_hyper);
            chain_V{i}=V;
            chain_V_L(i,:)=parametrization(V*v_hyper);
        else

        end
    end
    fprintf('\n')
    result.chain=chain;
    result.chain_V=chain_V;
    result.chain_V_L=chain_V_L;
    result.acceptance=accep;
end

%%
% function L_new=proposal(data,L_current)
%   L_new=mvnrnd(L_current,data.Sigma);
% end
function [L_new,Sigma,Sigma_new]=proposal(data,L_current,V,v)
    global MCMC_scale
    D_current=invparametrization(L_current,data.parameterization);
    d=size(D_current,1);  
    if nargin==2
        V=eye(d);
        v=d;
        prior_d=2*diag(d:-1:1);
        prior_offd=tril(ones(d,d),-1);
        prior_precision=prior_d+prior_offd';
        prior_precision=prior_precision(:);
    else
        prior_d=2*diag(v:-1:1);
        prior_offd=tril(ones(d,d),-1);
        prior_precision=prior_d+prior_offd';
        prior_precision=prior_precision*chol(V);
        prior_precision=prior_precision(:);
    end
    if d==2
    ind=[1,3,4];
    elseif d==3
    ind=[1,4,5,7,8,9];
    end
    prior_precision=prior_precision(ind);
    p=size(L_current,2);
    [y,var]=Emulator_Predict(data,D_current);
    n=size(var,1);
    Grad=gradient_tensor(data,L_current);
    Sigma=invandlogdet(Grad*((data.sigmasq*eye(n)+var)\Grad')+diag(1./prior_precision))*MCMC_scale;
    L_new=mvnrnd(L_current,Sigma);
    Grad_new=gradient_tensor(data,L_new);
    Sigma_new=invandlogdet(Grad_new*((data.sigmasq*eye(n)+var)\Grad_new')+diag(1./prior_precision))*MCMC_scale;
end

function y=gradient_tensor(data,L_current)
if data.iso_kernel==1
    D=data.para_emulator.D_para;
    q=size(D,1);
    C2=(sum((L_current-D).^2,2)).^(1/2);%distance_matrix(D,D_para,data.distance);
    eta=data.para_emulator.thetaopt(end);
    dR3=(-6/(eta^2))*exp(-sqrt(6)*C2./eta).*(ones(q,1)*L_current-D);
    y=dR3'*data.para_emulator.ResR2I';
else
    e=0.001;
    D_current=invparametrization(L_current,data.parameterization);
    y=Emulator_Predict(data,D_current);
    k=length(L_current);
    for i=1:k
        L_current1=L_current;
        L_current1(i)=L_current1(i)+e;
        D_current1=invparametrization(L_current1,data.parameterization);
        y1(i,:)=Emulator_Predict(data,D_current1);
    end
    y=(y1-ones(k,1)*y)/e;
end
end

function y=chol_bound(D)
 ee=eig(D);
 if max(ee)>=10
     y=1;
 else
     y=0;
 end
end