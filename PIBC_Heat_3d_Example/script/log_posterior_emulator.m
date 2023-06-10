function [y,var_emu]=log_posterior_emulator(data,D,V,v_n)
[y_s,var_emu]=Emulator_Predict(data,D);
n=size(var_emu,2);
Sigma=data.sigmasq*eye(n)+var_emu;
%Sigma=data.sigmasq*eye(n);
Rss=(y_s-data.obs)*(Sigma\(y_s-data.obs)');
prior=Prior_Density_PDT(D,V,v_n);
if var_emu(1,1)>0
    [s,v,d]=svd(var_emu/var_emu(1,1));
else
    v=zeros(n,n);
end
logdetsigma=sum(log(var_emu(1,1)*diag(v)+data.sigmasq));
y=-logdetsigma/2-Rss/2+log(prior);
end