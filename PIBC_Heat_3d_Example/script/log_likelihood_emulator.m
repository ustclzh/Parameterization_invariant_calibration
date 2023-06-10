function [y,var_emu]=log_likelihood_emulator(data,D)
[y_s,var_emu]=Emulator_Predict(data,D);
n=size(var_emu,2);
Sigma=data.sigmasq*eye(n)+var_emu;
%Sigma=data.sigmasq*eye(n);
Rss=(y_s-data.obs)*(Sigma\(y_s-data.obs)');

if var_emu(1,1)>0
    [s,v,d]=svd(var_emu/var_emu(1,1));
else
    v=zeros(n,n);
end
logdetsigma=sum(log(var_emu(1,1)*diag(v)+data.sigmasq));
y=-logdetsigma/2-Rss/2;
end