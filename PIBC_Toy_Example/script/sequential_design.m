function [D,x,hessian,data]=sequential_design(data,lik)
% if lik is assigned a value, then run the WLV design criterion
q_1=data.q_1;
dist=data.parameterization;
range=[0.01,-1,0.01,-1,-1,0.01;1,1,1,1,1,1]*10^(1/2);
dim=size(data.design_matrix,2);
range=range(:,1:dim*(dim+1)/2);
h=waitbar(0,'start');
for i=1:q_1
    str=[num2str(i),'of' num2str(q_1)];
    waitbar(i/q_1,h,str);
    global data_design
    data_design=data;
    options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000);
    options = optimoptions('particleswarm','Display','off','MaxIterations',3000,'HybridFcn',@patternsearch);
    if nargin==2
        [x01,fval0,exitflag0]= particleswarm(@WLV,dim*(dim+1)/2,range(1,:),range(2,:),options);
        [x2,fval2,exitflag2] = patternsearch(@WLV,[x01],[],[],[],[],range(1,:),range(2,:),[],options_patternsearch);
    else
        [x01,fval0,exitflag0]= particleswarm(@WPV,dim*(dim+1)/2,range(1,:),range(2,:),options);
        [x2,fval2,exitflag2] = patternsearch(@WPV,[x01],[],[],[],[],range(1,:),range(2,:),[],options_patternsearch);

    end

    D_add=invparametrization(x2,dist);
    dist_mat=pdist2(x2,data.design_para);
    if min(dist_mat)<0.05
        data.q_1=i-1;
        disp('early terminal for design generation')
        break
    end
    q=size(data.design_matrix,3);
    data.design_matrix(:,:,q+1)=D_add;
    data.design_para(q+1,:)=x2;
    data.Y_simulator(:,q+1)=simulator(data,D_add)';
    data=train_emulator(data,1);
end
close(h)
if nargin==2
[D,x,hessian]=map_optimization(data,lik);
else
[D,x,hessian]=map_optimization(data);
end
end


function y=WPV(L)
global data_design
%dist_to_design=min(pdist2(L,data_design.design_para));
D=invparametrization(L,data_design.parameterization);
[log_post,var_emu]=log_posterior_emulator(data_design,D);
y=-log_post-log(var_emu(1,1));%-log(dist_to_design);
end

function y=WLV(L)
global data_design
D=invparametrization(L,data_design.parameterization);
ee=eig(D);
if min(ee)<0.1 ||max(ee)>10
    y=inf;
    return
end
[log_post,var_emu]=log_likelihood_emulator(data_design,D);
y=-log_post-log(var_emu(1,1));
end
