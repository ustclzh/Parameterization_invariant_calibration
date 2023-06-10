function RESULT=Posterior_Infrerence(data_initial, D_true, strategy, do_mcmc, design_criterion)
% posterio inference
% strategy: 1:emulator, 2: simulator, 
% 3: emulator_likelihood, 4: simulator_likelihoos
% design_criterion: 1. for wlv, 2 for wpv
V = data_initial.V;
v = data_initial.v;

% optimization based on emulator+initial desin or simulator
RESULT.data_initial=data_initial;
[D_ini,x_ini,hessian_ini]=map_optimization(data_initial,strategy, V,v);
RESULT.err_ini=distance_matrix(D_ini,D_true,data_initial.parameterization);
RESULT.D_ini=D_ini;
RESULT.L_ini=x_ini;
RESULT.hessian_ini=hessian_ini;
RESULT.data_initial = data_initial;
if do_mcmc == 1
    RESULT.mcmc = mcmc_res(data_initial,D_ini,D_true);
end

% optimization based on emulator+initial desin+follow up design
if design_criterion >= 1 && strategy == 1
    data_seq=Design_Sequential(data_initial,V,v, design_criterion);
    [D_seq,x_seq,hessian_seq]=map_optimization(data_seq,strategy, V,v);
    RESULT.D_seq=D_seq;
    RESULT.L_seq=x_seq;
    RESULT.hessian_seq=hessian_seq;
    RESULT.data_seq=data_seq;
    RESULT.err_seq=distance_matrix(D_seq,D_true,data_initial.parameterization);
    if do_mcmc == 1
        RESULT.mcmc_seq = mcmc_res(data_seq,D_seq,D_true);
    end
end
%RESULT.bo1 = bayesian_optimization(data_initial,3, 30);
%RESULT.bo2 = bayesian_optimization_2(data_initial,3, 30);
RESULT.D_true=D_true;
end


function RESULT = mcmc_res(data_seq,D_seq,D_true)
    start_point = parametrization(D_seq,data_seq.parameterization);
    mcmc_seq = MCMC_tensor_emulator(data_seq,start_point,D_seq);
    D_seq_MCMC=invparametrization(MAP_est(mcmc_seq.chain(1000:end,:)));
    RESULT.D_MCMC=D_seq_MCMC;
    mcmc_seq.err=distance_matrix(D_seq_MCMC,D_true,data_seq.parameterization);
    RESULT.chain=mcmc_seq;
end


function y=MAP_est(sample)
    %return the MAP estimation for parameter, each column stands for a group of sample for one parameter, return the density mode of each
    %column.
    n=size(sample,2);
    for i=1:n
    data=sample(:,i);
    [a,b]=ksdensity(data);
    a_max=find(a==max(a));
    y(i)=b(a_max(1));
    end
end



