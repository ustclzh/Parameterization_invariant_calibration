function [y_pred,PosVar]=Emulator_Predict(data,D)
    para_emulator_all=data.para_emulator;
    y_pred = [];
    PosVar= [];
    iso_kernel=data.iso_kernel;
    if size(D,1) > 1 %if design is based on matrix, convert to vector parameterization
        D=parametrization(D,data.parameterization);
    end
    if size(D,1) == 1 %if design is based on matrix, convert to vector parameterization
        D_mat=invparametrization(D,data.parameterization);
    else
        D_mat = D;
    end 
    for emu_predict = 1 : length(data.data_structure)
        para_emulator=para_emulator_all{emu_predict};
        if iso_kernel==0 % product kernel
            D = parametrization(D_mat,data.distance_gp);
            D_in=para_emulator.D_in;
            D_para=para_emulator.D_para;
            cor_func=para_emulator.cor_func;
            thetaopt=para_emulator.thetaopt;
            alpha=para_emulator.alpha;
            sigma2=para_emulator.sigma2;
            R1=para_emulator.R1;
            ResR2I=para_emulator.ResR2I;
            R2I=para_emulator.R2I;
            [n1, d1]=size(D_in);
            R3=CompCorr(D,D_para,thetaopt((d1+1):end),cor_func)';
            dev=ResR2I*(R3');
            Ypred=alpha+dev;
            y_pred=[y_pred, Ypred'];%row vector
            PosVar = matrix_connection(PosVar, sigma2*R1*max(1+(1e-6)-R3*R2I*R3',0));%
            %PosVar = sigma2*R1*max(1+(1e-6)-R3*R2I*R3'+(1-M2TRI*R3')^2*MTRIMI,0);   
        elseif iso_kernel==1 % iso kernel
            D_in=para_emulator.D_in;
            D_para=para_emulator.D_para;
            D_matrix = para_emulator.D_matrix;
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
            n2 = size(D_para,1);
            %C=(sum((D-D_para).^2,2)).^(1/2);%distance_matrix(D,D_para,data.distance);
            C = zeros(n2,1);
            C(:,1) = distance_matrix(D_mat,D_matrix,data.distance_gp);
            R3=CompCorr_iso(C,thetaopt((d1+1):end),cor_func)';
            dev=ResR2I*R3';
            Ypred=alpha+dev;
            y_pred=[y_pred, Ypred'];%row vector
            PosVar=matrix_connection(PosVar, sigma2*R1*max(1+(1e-6)-R3*R2I*R3',0)); %
            %PosVar= sigma2*R1*max(1+(1e-6)-R3*R2I*R3'+(1-M2TRI*R3')^2*MTRIMI,0);   
        end
    end
end

function y = matrix_connection(A,B)
[n1, m1] = size(A);
[n2, m2] = size(B);
y = [A,zeros(n1,m2);zeros(n2,m1),B];
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
function Corr=CompCorr_iso(C2,eta,cor_func)
    if(cor_func==0)
        Corr=(eta.^((C2).^2));
    elseif(cor_func==1)
        rho1=sqrt(6)*C2./eta;
        Corr=((exp(-rho1)).*(rho1+1))+10^(-6)*(C2==0);
        Corr = prod(Corr,2);
    elseif(cor_func==2)
        rho1=2*sqrt(2.5)*C2./eta;
        Corr=((exp(-rho1)/3).*(rho1.^2+3*rho1+3)); 
    end
end