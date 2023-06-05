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
        case 3 % log-cholesky
            dim=floor((2*length(L))^(1/2));
            L = L -diag(diag(L)) + diag(exp(diag(L)));
            temp=zeros(dim,dim);
            ind=triu(ones(dim,dim));
            ind=find(ind>0);
            temp(ind)=L;
            D=temp'*temp;
        case 4 % spherical decomposition 
            dim=floor((2*length(L))^(1/2));
            ind=triu(ones(dim,dim));
            ind=find(ind>0);
            L_matrix = zeros(dim,dim);
            L_matrix(ind)=L;
            L = L_matrix;
            temp = zeros(dim,dim);
            for i =1:dim 
                temp (i,1) = L(i,1) * cos(L(i,2));
                for j=2:i
                    temp(i,j) = temp(i,j-1) / tan(L(i,j));
                end
            end
            index = 0:(dim^2-1);
            L = temp((rem(index,dim)<=floor(index/dim)));  
            D = L*L';
        case 5 % nothing
        

    end
end