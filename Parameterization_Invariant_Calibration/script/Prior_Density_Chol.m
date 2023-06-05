function y=Prior_Density_Chol(D,V,v)
% v = 1 proposed method 
% v = 0 noninformative
% return the density of wishart distribution given parameter V and v
%d=size(D,2);
%y=det(D)^((v-d-1)/2)*exp(-trace(V\D)/2);
if size(D,1) == 1 % if input a vector, then tranfer to matrix
    D = invparametrization(D);
end

if v == -1 % non-informative for cholesky
    prior_type = 1;
elseif v == 0 % non informative for matrix
    prior_type = 2;
elseif v == 1 % The matrix F prior
    prior_type = 3; 
else % whishart with specified hyperparameter
    prior_type = 4;
end
if min(eig(D))<=0
    y = 0;
    return
end
if prior_type == 1
    y = 1;
    return
elseif prior_type == 3
    d = size(D,1);
    [L,flag]=chol(D);
    y=det(D+2*eye(d))^(-(d+1)) * prod(diag(L)'.^(d:-1:1));
    return  
elseif prior_type == 2
    d = size(D,1);
    [L,flag]=chol(D);
    y=prod(diag(L)'.^(d:-1:1));
    return
else % whishart with specified hyperparameter
    d=size(D,2);
    [L,flag]=chol(V);
    [L,flag1]=chol((L'\D)/L);
    L1=diag(L)';
    L2=L-diag(L1);
    n_d=v+1-(1:d);
    y = prod(chi2pdf(L1.^2,n_d)).*prod(normpdf(L2(:)));
end

end