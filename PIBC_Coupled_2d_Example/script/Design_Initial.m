function data=Design_Initial(data,design_type,prior_type,V,v)
d=data.d;
if nargin == 1
    prior_type = 1;
    design_type = 1;
end
if nargin == 2
    prior_type = design_type;
end

data.design_type= design_type;
data.prior_type = prior_type;
if nargin >= 4
    data.V = V;
    data.v = v;
else
    switch prior_type % choose prior
        case 1 % non-informative for PDM
            data.V = eye(d);
            data.v = 0;
        case 2 % matrix-F
            data.V = eye(d);
            data.v = 1;
        case 3 % proposed prior, non-informative for chol
            data.V = eye(d);
            data.v = -1;
        case 4 % Wishart prior with specific V and v
            if nargin == 2
                data.V = eye(d)/(d+1);
                data.v = d+1;
            else
                data.V = V;
                data.v = v;
            end
    end
end


dist=data.parameterization;

% generating initial design
q_0=data.q_0;
range=[data.eig_L,data.eig_U];
tic
switch design_type %generate initial deisgn candidate
    case 1 % non informative
        b=rand_unif_PDT(10000*q_0,d,data.eig_U,data.distance_gp);
    case 2 % matrix F  
        b=rand_prop_PDT(10000*q_0,d,d*data.eig_U);
    case 3 % non-informative 
        b=rand_unif_chol(10000*q_0,d,d*data.eig_U);
    case 4 % wishart distribution
        b=rand_wish_PDT(10000*q_0,d,d*data.eig_U,V,v);
end
toc
% true value of theta is generated from the same as the prior type
                            switch prior_type %generate prior sample
                                case 1 % non informative
                                    temp=rand_unif_PDT(100000,d,data.eig_U,data.distance_gp);
                                case 2 % matrix F  
                                    temp=rand_prop_PDT(100000,d,d*data.eig_U);
                                case 3 % non-informative 
                                    temp=rand_unif_chol(100000,d,d*data.eig_U);
                                case 4 % wishart distribution
                                    temp=rand_wish_PDT(100000,d,d*data.eig_U,V,v);
                            end
                            
                            data.THETA_true=zeros(100000,d*(d+1)/2);
                            n=1;
                            for i=1:100000
                                D_0 = invparametrization(temp(i,:),dist);
                                [L,flag]=chol(D_0);
                                L = L';
                                if flag>0
                                    D_0
                                end
                                if max(eig(D_0))<=data.eig_U
                                    data.THETA_true(n,:)=temp(i,:);
                                    n=n+1;
                                end
                            end
                            data.THETA_true = data.THETA_true(1:n-1,:);
                            if n>1000
                                data.THETA_true = data.THETA_true(1:1000,:);
                            end



if d==2
    b(:,1)=abs(b(:,1));
    b(:,3)=abs(b(:,3));
elseif d==3
    b(:,1)=abs(b(:,1))+0.01;
    b(:,3)=abs(b(:,3))+0.01;
    b(:,6)=abs(b(:,6))+0.01;
end
c=zeros(10000*q_0,d*(d+1)/2);
n=1;
for i=1:10000*q_0
    D_0 = invparametrization(b(i,:),dist);
    [L,flag]=chol(D_0);
    L = L';
    if flag>0
        D_0
    end
    if max(eig(D_0))<=data.eig_U
        c(n,:)=b(i,:);
        n=n+1;
    end
end

c=c(1:n-1,:);
if n>1000*q_0 
    c = c(1:1000*q_0,:);
end
n
tic
%c = supportpoints(100*q_0,c,1000);
toc

if data.design_algo==1
    design=supportpoints(q_0,c,300);
elseif data.design_algo==2
    design=generate_design_finite_region(c,q_0,data); 
end
Design_matrix=zeros(d,d,q_0);
for i=1:q_0
   Design_matrix(:,:,i)=invparametrization(design(i,:),dist);
end
design_para=design(1:q_0,:);
data.design_para=design_para;
data.design_matrix=Design_matrix;
end

%%
function Y = rand_prop_PDT(n,d,r) %r upper bound of egenvalue
Y=zeros(n,d*(d+1)/2);
i = 1;
while i<=n
    V = iwishrnd(eye(d)*(2*d+2),d+1);
    X=wishrnd(V/(d+1),d+1);
    if max(eig(X))>r
        continue;
    end
    Y(i,:)=parametrization(X,1);
    i=i+1;
end
end


function Y = rand_wish_PDT(n,d,r, V, v) %r upper bound of egenvalue
Y=zeros(n,d*(d+1)/2);
i = 1;
while i<=n
    X=wishrnd(V,v);
    if max(eig(X))>r
        continue;
    end
    Y(i,:)=parametrization(X,1);
    i=i+1;
end
end

function Y = rand_unif_chol(n,d,r)%r upper bound of egenvalue
Y=zeros(n,d*(d+1)/2);
i = 1;
while i <=n
    seed = rand(1,d*(d+1)/2);
    for j =1:d*(d+1)/2
        if j == ceil(j^(1/2))*ceil(j^(1/2)+1)/2
            L(j) = r^(1/2)*seed(j);
        else
            L(j) = 2*r^(1/2)*seed(j)-r^(1/2);
        end
    end
    X = invparametrization(L,1);
    if max(eig(X))>r
        continue;
    end
    Y(i,:)=L;
    i= i+1;
end
end

function Y = rand_unif_PDT(n,d,r, distance_gp)%r upper bound of egenvalue
Y=zeros(n,d*(d+1)/2);
design = haltonset(d*(d+1)/2);
SEED = net(design,n+1);
SEED = SEED(2:end,:);
for i = 1:n
    seed = SEED(i,:);
    for j =1:d*(d+1)/2
        if distance_gp == 5
            if j == ceil(j^(1/2))*ceil(j^(1/2)+1)/2 % diagonal
                L(j) = r^(1/2)*seed(j)^(1/(d-floor((2*j)^(1/2))+1));
            else % off diagonal 
                L(j) = 2*r^(1/2)*seed(j)-r^(1/2);
            end
        elseif distance_gp == 1 %chol dist
            if j == ceil(j^(1/2))*ceil(j^(1/2)+1)/2 % diagonal
                L(j) = r^(1/2)*seed(j);
            else % off diagonal 
                L(j) = 2*r^(1/2)*seed(j)-r^(1/2);
            end
        end
    end
%    X = invparametrization(L,1);
%     if max(eig(X))>r
%         continue;
%     end
    Y(i,:)=L;
end
end

% function X=givens_rotation(d,i,j,theta)
%     X=eye(d);
%     X(i,j)=sin(theta);
%     X(i,i)=cos(theta);
%     X(j,i)=-sin(theta);
%     X(j,j)=cos(theta);
% 
% end

%%
function y=supportpoints(n,Sample,maxIte)
    % check
    % warning: chack Sample doesn't contain any repeatation
    % support points
    % n: desired sample size
    % Sample: candidate sample, i.e. MC sample of a target distribution
    % maxIte: max iteration, default 200;
    % if Sample consists of matrix, then return a design consists matrix.
    if length(size(Sample))>2
        res=1;
        Sample=spread(Sample);
    else
        res=0;
    end
    Sample=remove_rep(Sample);
    D_new=Sample(1:n,:);
    Sample(1:n,:)=[];
    p=size(Sample,2);
    D_0=zeros(n,p);
    if nargin==2
        maxIte=200;
    end
    ite=0;
    while sum(abs(D_0-D_new))>0.001
        
        ite=ite+1;
        D_0=D_new;
       for i=1:n
           D_new(i,:)=Mi(D_0,i,Sample);
       end
       if ite>=maxIte
           break;
       end
    end
    if res==1
        y=restore(D_new);
    else
        y=D_new;
    end
%     if checkingPDM(y)==1 
%     else
%         disp('wrong');
%     end
end

function y=Mi(D,i,sample)
    x=D(i,:);
    N=size(sample,1);
    p=length(x);
    n=size(D,1);
    d_2=sum((ones(N,1)*x-sample).^2,2).^(1/2);
    q=sum(d_2.^(-1));
    D(i,:)=[];
    d_x=sum((ones(n-1,1)*x-D).^2,2).^(1/2);
    D_w=(x-D)./(d_x*ones(1,p));
    x_new1=N*sum(D_w,1)/n;
    x_new2=sum(sample./(d_2*ones(1,p)),1);
    y=(x_new1+x_new2)/q;
end

function y=spread(sample)
    [p,q,n]=size(sample);
    y=zeros(n,p*q);
    for i=1:n
        temp=sample(:,:,i);
        y(i,:)=temp(:)';
    end
end

function y=restore(sample)
    [n,sqp]=size(sample);
    p=sqp^(1/2);
    y=zeros(p,p,n);
    for i=1:n
        temp=sample(i,:);
        y(:,:,i)=reshape(temp,p,p);
    end
end

function Sample=remove_rep(Sample)
    n=size(Sample,1);
    X=zeros(1,n);
    for i=1:n-1
        for j=i+1:n
            if sum(abs(Sample(i,:)-Sample(j,:)))==0
                X(j)=1;
            end
        end
    if sum(X)>=1
        Sample(X,:)=[];
    end
    end
end

function y=checkingPDM(sample)
    n=size(sample,3);
    y=1;
    for i=1:n
        v=eig(sample(:,:,i));
        if min(diag(v))<0
            y=0;
        end
    end
end
%%
function D_para=generate_design_finite_region(sample,d_desire,data)
    % check 1
    % design generated from prior sample
    % sample: candiate sample, sample size 2000;
    % d_desire: desired design size
    % if data.distance_gp == 5
    %     n = size(sample,1);
    %     d = size(sample,2);
    %     d = floor(sqrt(d*2));
    %     sample_cand = zeros(n,d^2);
    %     for i = 1:n
    %         temp = invparametrization(sample(i,:),data.parameterization);
    %         sample_cand(i,:) = temp(:);
    %     end
    % else
        sample_cand = sample;
    %end





    if size(sample_cand,1)>1000
    sample_cand=sample_cand(end-1000+1:end,:);
    end
    %
    % d=size(sample,1);
    % sample_para=zeros(size(sample,3),d*(d+1)/2);
    % for i=1:size(sample,3)
    %     D=sample(:,:,i);
    %     sample_para(i,:)=parametrization(D,data.parameterization);
    % end
    Design_index=Minimax_Design_from_Sample(sample_cand,d_desire,0.6,1,2,2);
    D_para=sample(Design_index,:);
end

function NearMinimaxSolution=Minimax_Design_from_Sample(CandidateSet,n,fl,fu,reduce,DP)
    %ProcedureC(CandidateSet,n,fl,fu,reduce,DP)
    %This program implements a heuristic procedure that finds near-minimax designs.
    %CandidateSet should be a matlab variable that contains the candidate set.
    %n=desired design size
    %L=fl*DKS, U=fu*DKS, DKS=distance of n-point KS design 
    %If reduce=1, then reduced SCLP is solved. Otherwise, SCLP is solved.
    %All distances are computed to DP decimal places.
    global N_fr StandardizedCS Distmat_fr 
    format long g
    if(n<=2||(ismember(mod(n,2),[0 1])~=1))
        disp('error')
        return
    end
    [N_fr, ~]=size(CandidateSet);
    StandardizedCS=(CandidateSet-repmat(mean(CandidateSet),N_fr,1))./repmat(std(CandidateSet,1,1),N_fr,1);
    Distmat_fr=pdist2(StandardizedCS,StandardizedCS);
    Distmat_fr=round(Distmat_fr*10^DP)/10^DP;
    [MaxDistforEachCol, RowIndices]=max(Distmat_fr);
    [~, MaxDistColumnIndex]=max(MaxDistforEachCol);
    MaxDistPair=[MaxDistColumnIndex RowIndices(MaxDistColumnIndex)];
    [~, InitialDesignIndices]=SequentialHeuristicDesign2(MaxDistPair,n);
    InitialNonDesignIndices=setdiff(1:N_fr,InitialDesignIndices);
    InitialDist=max(min(Distmat_fr(InitialDesignIndices,InitialNonDesignIndices)));
    U=fu*InitialDist;
    L=fl*InitialDist ;
    Phip=unique(Distmat_fr);
    m=length(Phip);
    theta=zeros(m,1);
    for i=1:(m-1)
        theta(i)=(Phip(i)+Phip(i+1))/2;
    end    
        theta(m)=1.1*Phip(m);
    index=(L<=Phip)&(Phip<=U);
    Sval=sort(theta(index),'descend');
    LS=length(Sval);
    Candidatedist=sort(Phip(index),'descend');
    storexopt=zeros(LS,N_fr);
    storez=zeros(LS,1);
    for i=1:LS
    if(i==1)
        check=1;
    else
        check=designdist>Sval(i);
    end
    if(check)
        indexmat=Distmat_fr<=Sval(i);
        if(reduce==1)
        [A, remainingvar]=trimmatrix(indexmat,N_fr);
        A=logical(A);
        else
        A=logical(indexmat);
        end
        xopt=HGW(A);
        z=sum(xopt);
    end
    storez(i)=z;
        if(reduce==1)
        storexopt(i,remainingvar)=xopt;  
        else
        storexopt(i,:)=xopt;
        end
    if(check)
    designdist=max(min(Distmat_fr(logical(storexopt(i,:)),logical(storexopt(i,:)==0))));
    end    
    end
    for i=1:LS-1
        for j=i+1:LS
            if(storez(i)>storez(j))
                storez(i)=storez(j);
            end
        end 
    end
    MaxDesignSizeIndex=storez==max(storez);
    storez(MaxDesignSizeIndex)=[];
    Candidatedist(MaxDesignSizeIndex)=[];
    MinimaxDesignSizes=unique(storez);
    NoDesignSizes=length(MinimaxDesignSizes);
    MinimaxDistances=zeros(NoDesignSizes,1);
    for i=1:NoDesignSizes
        IndexMinimaxDistance=find(abs(storez-MinimaxDesignSizes(i))<10^-5, 1, 'last' );
        MinimaxDistances(i)=Candidatedist(IndexMinimaxDistance);
    end
    DesignSize=find(MinimaxDesignSizes>=n,1);
    DesignSize=MinimaxDesignSizes(DesignSize);
    if isempty(DesignSize)
        DesignSize=MinimaxDesignSizes(end);
    end
    IndexNearMinimaxSolution=find(storez==DesignSize, 1, 'last' );    
    NearMinimaxSolution=storexopt(IndexNearMinimaxSolution,:)>(1-10^-5);
end

function [Design, DesignIndices]=SequentialHeuristicDesign2(StartDesignIndices,n)
    global StandardizedCS
    n0=length(StartDesignIndices);
    DesignIndices=StartDesignIndices;
    subset=StandardizedCS(StartDesignIndices,:);
    for i=1:(n-n0)
        NewPointIndex=findpoint(DesignIndices);
        DesignIndices=[DesignIndices NewPointIndex];
        subset=[subset; StandardizedCS(NewPointIndex,:)];
    end 
    Design=subset;
end  

function findpoint=findpoint(DesignIndices)
    global Distmat_fr N_fr
    NonDesignIndices=setdiff(1:N_fr,DesignIndices);
    Dist=min(Distmat_fr(DesignIndices,NonDesignIndices));
    [~, index]=max(Dist);
    findpoint=NonDesignIndices(index);
end

function HGW=HGW(A)
    [m, N]=size(A);
    b=-ones(m,1);
    c=ones(1,N);
    options=optimset('LargeScale','on','Display','off');
    [xopt , ~ , ~] = linprog(c,-A,b,[],[],zeros(1,N),ones(1,N),[],options);
    v=max(sum(A,2));
    tentative=xopt>=(1/v);
    Xi=find(tentative);
    n1=sum(tentative);
    stop=0;
    while(stop==0)
    Anew=A(:,Xi);
    redun=sum(Anew,2);
    mr=zeros(n1,1);
    for i=1:n1
        mr(i)=min(redun(Anew(:,i)));
    end
    [value, jstar]=max(mr);
    if(value>1)
        tentative(Xi(jstar))=0;
        Xi(jstar)=[];
        n1=n1-1;
    else
        stop=1;
    end
    end
    HGW=tentative;
end

function [trimmedmatrix, remainingvar]=trimmatrix(indexmat,N)
    i=1;
    n=N;
    while(i<n)
        check=0;
        j=i+1;
        while((check==0)&&(j<=n))
            Diff=indexmat(i,:)-indexmat(j,:);
            if(all(Diff>=0)&&(sum(Diff)>10^-6))
                indexmat(i,:)=[];
                n=n-1;
                check=1;
            elseif(all(Diff<=0))
                indexmat(j,:)=[];
                n=n-1;      
            else
                j=j+1;
            end
        end
        if(check==0)
            i=i+1;
        end
    end
    
    storeindex=1:N;
    i=1;
    n=N;
    while(i<n)
        check=0;
        j=i+1;
        while((check==0)&&(j<=n))
            Diff=indexmat(:,j)-indexmat(:,i);
            if(all(Diff>=0)&&(sum(Diff)>10^-6))
                indexmat(:,i)=[];
                storeindex(i)=[];
                n=n-1;
                check=1;
            elseif(all(Diff<=0))
                indexmat(:,j)=[];
                storeindex(j)=[];
                n=n-1;                
            else
                j=j+1;
            end
        end
        if(check==0)
            i=i+1;
        end
    end
    trimmedmatrix=indexmat;
    remainingvar=storeindex;
end

