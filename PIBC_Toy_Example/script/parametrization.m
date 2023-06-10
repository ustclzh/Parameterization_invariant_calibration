function L=parametrization(D,parameterization)
if nargin==1
    parameterization=1;
end
switch parameterization
    case 1
        dim=size(D,1);
        temp=chol(D);
        temp=temp(:)';
        index=0:(dim^2-1);
        L=temp((rem(index,dim)<=floor(index/dim)));
    case 2
        dim=size(D,1);
        temp=logm(D);
        temp=temp(:)';
        index=0:(dim^2-1);
        L=temp((rem(index,dim)<=floor(index/dim)));
end
end