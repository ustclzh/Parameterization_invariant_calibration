function result=MCMC_tensor(data,para_MCMC,start_point,Sigma)
    N=para_MCMC.N;
    D=invparametrization(start_point,data.parameterization);
    if max(abs(log(eig(D))))>10
        start_point=data.design_matrix(:,:,1);
    end
    M=length(start_point);
    chain=zeros(N,M);
    L_current=start_point;
    k=size(L_current,2);
    D_current=invparametrization(L_current,data.parameterization);
    logpost_current=log_posterior_simulator(data,D_current);
    accep=0;
    Grad=gradient_tensor(data,L_current);
    n=size(Grad,2);
    if nargin==3 % Specify stating point
        data.Sigma=invandlogdet(Grad*((data.sigmasq*eye(n))\Grad')+eye(k));
    else
        data.Sigma=Sigma;
    end
    h=waitbar(0,'start');
    time_used=0;
    for i=1:N
        if nargin==3
            if rem(i,100)==0
                Grad=gradient_tensor(data,L_current);
                data.Sigma=invandlogdet(Grad*((data.sigmasq*eye(n))\Grad')+eye(k));
            end
        end
        tic
        L_new=proposal(data,L_current);
        D_new=invparametrization(L_new,data.parameterization);
        logpost_new=log_posterior_simulator(data,D_new);
        alpha=exp(logpost_new-logpost_current);
        
        u=rand(1);
        if alpha>u
            chain(i,:)=L_new;
            L_current=L_new;
            accep=accep+1;
            logpost_current=logpost_new;
        else
            chain(i,:)=L_current;
        end
        time_Z=toc;
        time_used=time_used+time_Z;
        if rem(i,100)==0
            str=[num2str(floor(N/100-i/100)) '00 steps, ' num2str(floor(time_used*(N-i)/i)), 's.',' P_a: ' num2str(floor(100*accep/i))];
            waitbar(i/N,h,str);
        end
    end
    close(h);
    result.chain=chain;
    result.acceptance=accep;
end
%%
function L_new=proposal(data,L_current)
    L_new=mvnrnd(L_current,2*data.Sigma);
end

function y=gradient_tensor(data,L_current)
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

