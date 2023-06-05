function [y_pred,PosVar]=emulator(data,D)
    para_emulator=data.para_emulator;
    iso_kernel=data.iso_kernel;
    if size(D,1)>1
        D=parametrization(D,data.parameterization);
    end
    if iso_kernel==0 % product kernel
        D_in=para_emulator.D_in;
        D_para=para_emulator.D_para;
        cor_func=para_emulator.cor_func;
        thetaopt=para_emulator.thetaopt;
        alpha=para_emulator.alpha;
        sigma2=para_emulator.sigma2;
        R1=para_emulator.R1;
        ResR2I=para_emulator.ResR2I;
        M2TRI=para_emulator.M2TRI;
        MTRIMI=para_emulator.MTRIMI;
        R2I=para_emulator.R2I;
        [n1, d1]=size(D_in);
        R3=CompCorr(D,D_para,thetaopt((d1+1):end),cor_func)';
        dev=ResR2I*(R3');
        Ypred=alpha+dev;
        y_pred=Ypred';%row vector
        PosVar = sigma2*R1*max(1+(1e-6)-R3*R2I*R3',0);%
        %PosVar = sigma2*R1*max(1+(1e-6)-R3*R2I*R3'+(1-M2TRI*R3')^2*MTRIMI,0);   
    elseif iso_kernel==1 % iso kernel
        D_in=para_emulator.D_in;
        D_para=para_emulator.D_para;
        cor_func=para_emulator.cor_func;
        thetaopt=para_emulator.thetaopt;
        alpha=para_emulator.alpha;
        sigma2=para_emulator.sigma2;
        M2TRI=para_emulator.M2TRI;
        MTRIMI=para_emulator.MTRIMI;
        R1=para_emulator.R1;
        ResR2I=para_emulator.ResR2I;
        R2I=para_emulator.R2I;
        [n1, d1]=size(D_in);
        C=(sum((D-D_para).^2,2)).^(1/2);%distance_matrix(D,D_para,data.distance);
        R3=CompCorr2(C',thetaopt((d1+1):end),cor_func);
        dev=ResR2I*R3';
        Ypred=alpha+dev;
        y_pred=Ypred';%row vector
        PosVar=sigma2*R1*max(1+(1e-6)-R3*R2I*R3',0); %
        %PosVar= sigma2*R1*max(1+(1e-6)-R3*R2I*R3'+(1-M2TRI*R3')^2*MTRIMI,0);   
    end
end


function y=CompCorr(x1,x2,eta,cor_func)
    if(cor_func==0)
        y=prod(eta.^((x1-x2).^2),2);
    elseif(cor_func==1)
        rho=x1-x2;
        rho1=sqrt(6)*abs(rho)./eta;
        y=prod((exp(-rho1)).*(rho1+1),2);
    elseif(cor_func==2)
        rho=x1-x2;
        rho1=2*sqrt(2.5)*abs(rho)./eta;
        y=prod((exp(-rho1)/3).*(rho1.^2+3*rho1+3),2); 
    end
end
function Corr=CompCorr2(C2,eta,cor_func)
    if(cor_func==0)
        Corr=(eta.^((C2).^2));
    elseif(cor_func==1)
        rho1=sqrt(6)*C2./eta;
        Corr=((exp(-rho1)).*(rho1+1))+10^(-6)*(C2==0);
    elseif(cor_func==2)
        rho1=2*sqrt(2.5)*C2./eta;
        Corr=((exp(-rho1)/3).*(rho1.^2+3*rho1+3)); 
    end
end