function data=Design_Sequential(data,V,v,design_criterion)
% if lik is assigned a value, then run the WLV design criterion
q_1=data.q_1;
data.map_trace = [];
data.err_trace = [];
data.y_err_trace = [];
dim=size(data.design_matrix,2);
range=design_region(data);

fprintf('\n Sequential Design ')
for i=1:q_1
    if  i> 1
    fprintf(repmat('\b',1,10))
    end
    fprintf('%3d of %3d',i,q_1)
    f_obj = @(L) WLV(data, L);
    strategy = 1;
    if design_criterion == 2
        f_obj = @(L) WPV(data, L, V, v);
    end
    options_patternsearch = optimoptions('patternsearch','Display','off','MaxIterations',3000);
    options = optimoptions('particleswarm','Display','off','MaxIterations',3000,'HybridFcn',@patternsearch);
    [x01,fval0,exitflag0]= particleswarm(f_obj,dim*(dim+1)/2,range(1,:),range(2,:),options);
    [x2,fval2,exitflag2] = patternsearch(f_obj,[x01],[],[],[],[],range(1,:),range(2,:),[],options_patternsearch);
    options_neal=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter', 5000, 'Display', 'off');
    [x3,fval3]=fminsearch(f_obj,x2,options_neal);

    L_true = parameterization(data.D_true);
    [x4,fval4]=fminsearch(f_obj,L_true,options_neal);
    if fval4 <=fval3
    x3 = x4;
    fval3 = fval4;
    end

    if fval3<fval2
        x2 = x3;
    end



    D_add=invparametrization(x2,data.parameterization);
    dist_mat=pdist2(x2,data.design_para);
    if min(dist_mat)<1e-8
        data.q_1=i-1;
        disp('early termination for design generation')
        break
    end
    q=size(data.design_matrix,3);
    data.design_matrix(:,:,q+1)=D_add;
    data.design_para(q+1,:)=x2;
    y_add = simulator(data,D_add);
    
    data.Y_simulator(:,q+1)=y_add';
    data=Emulator_Train(data,1,0);

    data.y_err_trace = [data.y_err_trace,mean((y_add-data.obs).^2)];
    if rem(i,data.d*(data.d+1)/2) == 0
        [D_ini,x_ini,hessian_ini]=map_optimization(data,strategy, V, v);
        data.map_trace = [data.map_trace;x_ini];
        data.err_trace = [data.err_trace, distance_matrix(D_ini,data.D_true,data.parameterization)];
    end


end
fprintf('\n')
end


function y=WPV(data, L, V,v)
D=invparametrization(L,data.parameterization);
e = eig(D);
if min(e)<data.eig_L || max(e)>data.eig_U
    y = 10^256;
    return;
end
[log_post,var_emu]=log_posterior_emulator(data,D,V,v);
y=-log_post-log(max(diag(var_emu)));
end

function y=WLV(data, L)
D=invparametrization(L,data.parameterization);
e = eig(D);
if min(e)<data.eig_L || max(e)>data.eig_U
    y = 10^256;
    return;
end
[log_post,var_emu]=log_likelihood_emulator(data,D);
y=-log_post-log(max(diag(var_emu)));
end