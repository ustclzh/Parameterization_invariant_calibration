function L=parametrization(D,parameterization)
if nargin==1
    parameterization=1;
end
switch parameterization
    case 1 %cholesky
        dim=size(D,1);
        [temp,flag]=chol(D);
        if flag > 0 
            L = zeros(1,dim*(dim+1)/2);
            return;
        end
        temp=temp(:)';
        index=0:(dim^2-1);
        L=temp((rem(index,dim)<=floor(index/dim)));
    case 2 % matrix logrithm
        dim=size(D,1);
        temp=logm(D);
        temp=temp(:)';
        index=0:(dim^2-1);
        L=temp((rem(index,dim)<=floor(index/dim)));
    case 3 % log-cholesky
        dim = size(D,1);
        [temp,flag] = chol(D);
        temp = temp -diag(diag(temp)) + diag(log(diag(temp)));
        temp = temp(:)';
        index = 0:(dim^2-1);
        L = temp((rem(index,dim)<=floor(index/dim)));  
    case 4 % spherical decomposition 
        dim = size(D,1);
        [temp,flag] = chol(D);
        L = zeros(dim,dim);
        for i =1:dim 
            L(i,1) = sum(temp(i,:).^2)^(1/2);
            for j=2:i
                L(i,j)=asin((L(i,1)^2 - sum(temp(i,1:j-1).^2)^(1/2)/temp(i,j)));
            end
        end
        temp = L(:)';
        index = 0:(dim^2-1);
        L = temp((rem(index,dim)<=floor(index/dim)));  
    case 5 %original 
        dim=size(D,1);
        temp=D(:)';
        index=0:(dim^2-1);
        L=temp((rem(index,dim)<=floor(index/dim)));

end
end