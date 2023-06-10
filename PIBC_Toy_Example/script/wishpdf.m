function y=wishpdf(D,V,v)
% return the density of wishart distribution given parameter V and v
%d=size(D,2);
%y=det(D)^((v-d-1)/2)*exp(-trace(V\D)/2);
if nargin ==1
    V=eye(size(D,1));
    v=size(D,1);
    if max(eig(D))>10 || min(eig(D))<0.1
        y=0;
        return
    end
    L=chol(D);
    y=prod(diag(L)'.*(v:-1:1));
    return
end
if max(eig(D))>10 || min(eig(D))<0.1
    y=0;
    return
end
d=size(D,2);
L=chol(V);
L=chol((L'\D)/L);
L1=diag(L)';
L2=L-diag(L1);
n_d=v+1-(1:d);
y=prod(chi2pdf(L1.^2,n_d)).*prod(normpdf(L2(:)))/(normpdf(0)^(d*(d+1)/2));
%y=1;
end