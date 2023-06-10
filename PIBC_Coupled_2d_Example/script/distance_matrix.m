function y=distance_matrix(A,B,dist)
% return the distance of two matrix A and B;
% A, B are two designs, return a distance matrix whose i,j element is the
% distance between A(:,:,i) and B(:,:,j)
% dist: the distance metric
% 1: cholesky distance;
% 2: matrix logarithm distance;
% 3: Geodetic distance ||log(eig(A^-1B))||;
% 4: log cholesky distance.
if nargin == 2
    dist = 5;
end

N1=size(A);
N2=size(B);
if length(N1)==2
    N1(3)=1;
    A(:,:,1)=A;
end
if length(N2)==2
    N2(3)=1;
    B(:,:,1)=B;
end
y=zeros(N1(3),N2(3));
switch dist    
    case 1% chol 
        for i=1:N1(3)
            for j=1:N2(3)
                [L_a,flag_a]=chol(A(:,:,i));
                [L_b,flag_b]=chol(B(:,:,j));
                temp=sum(sum((L_a-L_b).^2))^(1/2);
                y(i,j)=temp;
                %y(i,j)=dist_mat(A(:,:,i),B(:,:,j));
            end
        end
    case 2 %logm
        for i=1:N1(3)
            for j=1:N2(3)
                L_a=tril(logm(A(:,:,i)));
                L_b=tril(logm(B(:,:,j)));
                temp=sum(sum((L_a-L_b).^2))^(1/2);
                y(i,j)=temp;
            end
        end  
    case 3 % Geodetic distance
        for i=1:N1(3)
            for j=1:N2(3)
                a=eig(A(:,:,i)\B(:,:,j));
                y(i,j)=(sum((log(a)).^2))^(1/2);
            end
        end  
    case 4 % log chol
        for i=1:N1(3)
            for j=1:N2(3)
                [L_a,flag_a]=chol(A(:,:,i));
                [L_b,flag_b]=chol(B(:,:,j));
                L_a=tril(L_a,-1)+diag(log(diag(L_a)));
                L_b=tril(L_b,-1)+diag(log(diag(L_b)));
                temp=sum(sum((L_a-L_b).^2))^(1/2);
                y(i,j)=temp;
            end
        end
    case 5 % L2
         for i=1:N1(3)
            for j=1:N2(3)
                L_a=A(:,:,i);
                L_b=B(:,:,j);
                temp=sum(sum((L_a-L_b).^2))^(1/2);
                y(i,j)=temp;
            end
        end
end
end