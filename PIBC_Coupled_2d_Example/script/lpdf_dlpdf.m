
function [lpdf,dlpdf]=lpdf_dlpdf(data,L,V,v_n)
parameterization = data.parameterization;
D = invparametrization(L,parameterization);
lpdf=lpdf_cal(data,D,V,v_n);
k = length(L);
e = 0.0001;
dlpdf = zeros(1,k);
for i =1:k
K = L;
K(i) = K(i)+e;
D = invparametrization(K,parameterization);
dlpdf(i)=lpdf_cal(data,D,V,v_n);
end
end


function lpdf=lpdf_cal(data,D,V,v_n)
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
lpdf=-logdetsigma/2-Rss/2+log(prior);
end