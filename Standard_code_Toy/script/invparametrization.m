function D=invparametrization(L,parameterization)
    if nargin==1
        parameterization=1;
    end
    switch parameterization
        case 1
            dim=floor((2*length(L))^(1/2));
            temp=zeros(dim,dim);
            ind=triu(ones(dim,dim));
            ind=find(ind>0);
            temp(ind)=L;
            D=temp'*temp;
        case 2
            dim=floor((2*length(L))^(1/2));
            temp=zeros(dim,dim);
            ind=triu(ones(dim,dim));
            ind=find(ind>0);
            temp(ind)=L;
            D=temp+temp';
            D=D-diag(diag(D)/2);
            D=expm(D);
    end
end