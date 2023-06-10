function data=Emulator_Train(data,indiator, update_hyper)
% check
% indicator: if nonempty, means to update emulator only, i.e. no need to add simulation data
k = length(data.data_structure); % number of PDE in the equation
%disp(['number of simulator:',num2str(k)])
if k>=2
    useage_ind = [];
    for i =1:k
        useage_ind = matrix_connection(useage_ind,ones(1,data.data_structure(i)));
    end
else
    useage_ind = ones(1,data.data_structure(1));
end

if nargin <=2
    update_hyper = 1;
end

data.para_emulator = cell(1,k);
if nargin==1 % need to get simulator output. 
    data=outputfromdesign(data);
end
for emu_train = 1:k
    D_in=data.D_in(useage_ind(emu_train,:)>0,:);
    D=data.design_matrix;
    D_para=data.design_para;
    iso_kernel=data.iso_kernel;
    if iso_kernel==1 % isotropic kernel
        %D_para=data.design_matrix;
        cor_func=data.cor_func;
        No=cor_func;
        Y_temp=data.Y_simulator(useage_ind(emu_train,:)>0,:);
        [n1, d1]=size(D_in);
        n2=size(D_para,1);
        n=n1*n2;
        d=d1+1;
        d2 = 1;
        if(sum(abs(size(Y_temp)-[n1 n2]))>10^-6)
            display('error')
            return    
        end
        scale=[range(D_in,1),10*ones(1,1)];
        C1=cell(1,d1);
        for i=1:d1
            C1{i}=abs(repmat(D_in(:,i),1,n1)-repmat(D_in(:,i)',n1,1));
        end
        C2{1}=distance_matrix(D,D,data.distance_gp);
        options=optimoptions(@fminunc,'MaxIter',10^3,'TolX',0.01,'TolFun',0.01,'MaxFunEvals',10^5,'Display','off','Algorithm','quasi-newton','ObjectiveLimit',-10^250);
        Obj = @(par) Lkh(par, Y_temp, n, d, No, n1, n2, d1, d2, scale,C1,C2);
        if(No==0)
            [paropt1, fval1, exitflag1]=fminunc(Obj,[norminv(0.25.^(1./(scale.^2)))],options); 
            [paropt2, fval2, exitflag2]=fminunc(Obj,[norminv(0.75.^(1./(scale.^2)))],options);
        elseif(No==1)
            [paropt1, fval1, exitflag1]=fminunc(Obj,[log(0.909700041540068.*scale)],options);  
            [paropt2, fval2, exitflag2]=fminunc(Obj,[log(2.54815755509488.*scale)],options);       
        elseif(No==2)
            [paropt1, fval1, exitflag1]=fminunc(Obj,log(0.888569182612929.*scale),options); 
            [paropt2, fval2, exitflag2]=fminunc(Obj,log(2.24172067901287.*scale),options);    
        end
        if((exitflag1<=0)||(exitflag2<=0))
            display('Likelihood optimization failure: algorithm did not converge.')   
        end
        if(fval1<fval2)
            if(No==0)
                thetaopt=normcdf(paropt1);
            elseif(No==1)
                thetaopt=exp(paropt1);        
            elseif(No==2)
                thetaopt=exp(paropt1);        
            end
        else
            if(No==0)
                thetaopt=normcdf(paropt2);
            elseif(No==1)
                thetaopt=exp(paropt2);  
            elseif(No==2)
                thetaopt=exp(paropt2);          
            end
        end
        %display('Optimal values of GP parameters')
        %disp(thetaopt)
        R1=CompCorr1(thetaopt(1:d1),No, C1, d1);
        %R1=eye(n1);
        R2=CompCorr2(thetaopt((d1+1):d),No, C2, d2);
        R1I=invandlogdet(R1);
        R2I=invandlogdet(R2);
        M2TRI=ones(1,n2)*R2I;
        SSR2I=M2TRI*ones(n2,1);
        %alpha=mean(Y_temp,2);%
        alpha=Y_temp*(M2TRI')/SSR2I;
        MTRIMI=R1/SSR2I;
        Res=Y_temp-alpha*ones(1,n2);
        A=Res'*R1I;
        ResR2I=Res*R2I;
        %sigma2=sum(sum(A'.*ResR2I))/n;%
        sigma2=sum(sum(A'.*ResR2I))/(n-n2);
        para_emulator.D_in=D_in;
        para_emulator.D_para=D_para;
        para_emulator.D_matrix = D;
        para_emulator.cor_func=cor_func;
        para_emulator.thetaopt=thetaopt;
        para_emulator.alpha=alpha;
        para_emulator.sigma2=sigma2;
        para_emulator.R1=R1;
        para_emulator.R2=R2;
        para_emulator.ResR2I=ResR2I;
        para_emulator.R2I=R2I;
        para_emulator.MTRIMI=MTRIMI;
        para_emulator.M2TRI=M2TRI;
        para_emulator.Y=Y_temp;
        data.para_emulator{emu_train}=para_emulator;
    elseif iso_kernel==0 % product kernel
        D_para=data.design_para;
        cor_func=data.cor_func;
        No=cor_func;
        Y_temp=data.Y_simulator(useage_ind(emu_train,:)>0,:);
        [n1, d1]=size(D_in);
        [n2, d2]=size(D_para);
        for i =1:n2
            D_para(i,:) = parametrization(data.design_matrix(:,:,i),data.distance_gp);
        end
        data.design_para = D_para;
        
        n=n1*n2;
        d=d1+d2;
        if(sum(abs(size(Y_temp)-[n1 n2]))>10^-6)
            display('error')
            return    
        end
        scale=[range(D_in,1) range(D_para,1)];
    
        C1=cell(1,d1);
        for i=1:d1
            C1{i}=abs(repmat(D_in(:,i),1,n1)-repmat(D_in(:,i)',n1,1));
        end
        C2=cell(1,d2);
        for i=1:d2
            C2{i}=abs(repmat(D_para(:,i),1,n2)-repmat(D_para(:,i)',n2,1));
        end
        options=optimoptions(@fminunc,'MaxIter',10^3,'TolX',0.01,'TolFun',0.01,'MaxFunEvals',10^5,'Display','off','Algorithm','quasi-newton','ObjectiveLimit',-10^250);
        Obj = @(par) Lkh(par, Y_temp, n, d, No, n1, n2, d1, d2, scale,C1,C2);
        if(No==0)
            [paropt1, fval1, exitflag1]=fminunc(Obj,[norminv(0.25.^(1./(scale.^2)))],options); 
            [paropt2, fval2, exitflag2]=fminunc(Obj,[norminv(0.75.^(1./(scale.^2)))],options);
        elseif(No==1)
            [paropt1, fval1, exitflag1]=fminunc(Obj,[log(0.909700041540068.*scale)],options);
            [paropt2, fval2, exitflag2]=fminunc(Obj,[log(2.54815755509488.*scale)],options);  
        elseif(No==2)
            [paropt1, fval1, exitflag1]=fminunc(Obj,log(0.888569182612929.*scale),options); 
            [paropt2, fval2, exitflag2]=fminunc(Obj,log(2.24172067901287.*scale),options);    
        end
        if((exitflag1<=0)||(exitflag2<=0))
            display('Likelihood optimization failure: algorithm did not converge.')   
        end
        if(fval1<fval2)
            if(No==0)
                thetaopt=normcdf(paropt1);
            elseif(No==1)
                thetaopt=exp(paropt1);        
            elseif(No==2)
                thetaopt=exp(paropt1);        
            end
        else
            if(No==0)
                thetaopt=normcdf(paropt2);
            elseif(No==1)
                thetaopt=exp(paropt2);  
            elseif(No==2)
                thetaopt=exp(paropt2);          
            end
        end
        
%       display('Optimal values of GP parameters')
%       disp(thetaopt)
        R1=CompCorr1(thetaopt(1:d1),No, C1, d1);
        R2=CompCorr2(thetaopt((d1+1):d),No, C2, d2);
        R1I=invandlogdet(R1);
        R2I=invandlogdet(R2);
        M2TRI=ones(1,n2)*R2I;
        SSR2I=M2TRI*ones(n2,1);
        %alpha=mean(Y_temp,2);%
        alpha=Y_temp*(M2TRI')/SSR2I;
        MTRIMI=R1/SSR2I;
        Res=Y_temp-alpha*ones(1,n2);
        A=Res'*R1I;
        ResR2I=Res*R2I;
        %sigma2=sum(sum(A'.*ResR2I))/n;%/(n-n2)
        sigma2=sum(sum(A'.*ResR2I))/(n-n2);
        para_emulator.D_in=D_in;
        para_emulator.D_para=D_para;
        para_emulator.cor_func=cor_func;
        para_emulator.thetaopt=thetaopt;
        para_emulator.alpha=alpha;
        para_emulator.sigma2=sigma2;
        para_emulator.R1=R1;
        para_emulator.R2=R2;
        para_emulator.ResR2I=ResR2I;
        para_emulator.R2I=R2I;
        para_emulator.MTRIMI=MTRIMI;
        para_emulator.M2TRI=M2TRI;
        para_emulator.Y=Y_temp;
        data.para_emulator{emu_train}=para_emulator;
    end
end
end
%%
function y = matrix_connection(A,B)
    [n1, m1]= size(A);
    [n2, m2] = size(B);
    y = [A,zeros(n1,m2);zeros(n2,m1),B];
end

function Objective=Lkh(par, Y_temp, n, d, No, n1, n2, d1, d2, scale,C1,C2)
    if(sum(isnan(par))>0)
        Objective=Inf;
        return
    end
    if(No==0)
        r=normcdf(par);
        if(sum(r>(0.999.^(1./(scale.^2))))>0)
            Objective=Inf;
            return
        end
    elseif(No==1)
        r=exp(par);    
        if(sum(r>(53.9511207457682.*scale))>0)
            Objective=Inf;
            return        
        end      
    elseif(No==2)
        r=exp(par);
        if(sum(r>(40.7953930912638.*scale))>0)
            Objective=Inf;
            return        
        end
    end
    R1=CompCorr1(r(1:d1),No, C1, d1);
    %R1=eye(n1);
    R2=CompCorr2(r((d1+1):d),No, C2, d2);
    [R1I,LD1]=invandlogdet(R1);
    [R2I,LD2]=invandlogdet(R2);
    M2TRI=ones(1,n2)*R2I;
    SSR2I=M2TRI*ones(n2,1);
    LDR=(n2)*LD1+(n1)*LD2;
    %alpha=mean(Y_temp,2);
    alpha=Y_temp*(M2TRI')/SSR2I;
    Res=Y_temp-alpha*ones(1,n2);
    A=Res'*R1I;
    B=Res*R2I;
    %sigma2=sum(sum(A'.*B))/n;%
    sigma2=sum(sum(A'.*B))/(n-n2);
    Objective=n*log(max(sigma2,0))+LDR;
    if(Objective<(-10^250))
        Objective=-10^250;
    elseif(isnan(Objective)||(Objective>10^250))
        Objective=10^250;    
    end
end

function Corr=CompCorr1(r,No, C1, d1)
    if(No==0)
        Corr=1;
        for i=1:d1
            Corr=Corr.*(r(i).^((C1{i}).^2));
        end
    elseif(No==1)
	    Corr=1;
        for i=1:d1
            rho1=sqrt(6)*C1{i}./r(i);
            Corr=Corr.*((exp(-rho1)).*(rho1+1));
        end
    elseif(No==2)
        Corr=1;
        for i=1:d1
            rho1=2*sqrt(2.5)*C1{i}./r(i);
            Corr=Corr.*((exp(-rho1)/3).*(rho1.^2+3*rho1+3)); 
        end
    end
    Corr=Corr+10^(-6)*eye(size(Corr,1));
end

function Corr=CompCorr2(r,No, C2, d2)
    if(No==0)
        Corr=1;
        for i=1:d2
            Corr=Corr.*(r(i).^((C2{i}).^2));
        end
    elseif(No==1)
	    Corr=1;
        for i=1:d2
            rho1=sqrt(6)*C2{i}./r(i);
            Corr=Corr.*((exp(-rho1)).*(rho1+1));
        end
    elseif(No==2)
        Corr=1;
        for i=1:d2
            rho1=2*sqrt(2.5)*C2{i}./r(i);
            Corr=Corr.*((exp(-rho1)/3).*(rho1.^2+3*rho1+3)); 
        end
    end
    Corr=Corr+10^(-6)*eye(size(Corr,1));
end
    
function data=outputfromdesign(data)
    % design is a 3 dimensional variable, design(:,:,i) is the  i-th design
    % point.
    % output: each  row stands for simulation output of i-th design point.
    n=size(data.design_matrix,3);
    for i=1:n
       y(i,:)=simulator(data,data.design_matrix(:,:,i));
    end
    data.Y_simulator=y';
end