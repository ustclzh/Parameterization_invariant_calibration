function RESULT=run_emulator_based_infer_full(data_initial,data_initial1,D_true)
%%
    [D_ini,x_ini,hessian_ini]=map_optimization(data_initial);
    disp('MAP of initial design based Posterior and err')
    D_ini
    Err(1)=distance_matrix(D_ini,D_true,data_initial.parameterization)

    [D_ini_h,x_ini_h,hessian_ini_h]=map_optimization(data_initial,1);
    disp('MAP of initial design based Posterior and err') 
    D_ini_h
    Err(2)=distance_matrix(D_ini,D_true,data_initial.parameterization)
    
    para_MCMC.N = data_initial.N_MCMC;
    start_point = parametrization(D_ini,data_initial.parameterization);
    result_emu_initial = MCMC_tensor_emulator(data_initial,para_MCMC,start_point);
    result_emu_initial_h = MCMC_tensor_emulator_hierarchical(data_initial,para_MCMC,start_point);
    

    % optimization based on emulator+initial desin+follow up design
    [D_seq,x_seq,hessian_seq,data_seq]=sequential_design(data_initial);

    disp('MAP of sequential design based Posterior and err')
    D_seq 
    Err(3)=distance_matrix(D_seq,D_true,data_initial.parameterization)

    para_MCMC.N = data_initial.N_MCMC;
    start_point = parametrization(D_seq,data_initial.parameterization);
    result_emu_seq = MCMC_tensor_emulator(data_seq,para_MCMC,start_point);
    %result_emu_seq_h = MCMC_tensor_emulator_hierarchical(data_seq,para_MCMC,start_point);
    % optimization based on emulator+initial desin+follow up design
    [D_seq_lik,x_seq,hessian_seq_lik,data_seq_lik]=sequential_design(data_initial,1);

    disp('MAP of sequential design based Posterior and err')
    D_seq_lik 
    Err(4)=distance_matrix(D_seq_lik,D_true,data_initial.parameterization)

    para_MCMC.N = data_initial1.N_MCMC;
    start_point = parametrization(D_seq_lik,data_initial.parameterization);
    result_emu_seq_lik = MCMC_tensor_emulator_hierarchical(data_seq_lik,para_MCMC,start_point);


    % optimization based on emulator+support point+follow up design
    [D_ini_only,x_ini_only,hessian_ini_only]=map_optimization(data_initial1);
    disp('MAP of SP design based Posterior and err')
    D_ini_only 
    Err(5)=distance_matrix(D_ini_only,D_true,data_initial.parameterization)

    [D_ini_only_h,x_ini_only_h,hessian_ini_only_h]=map_optimization(data_initial1,1);
    disp('MAP of SP design based Posterior and err')
    D_ini_only 
    Err(6)=distance_matrix(D_ini_only,D_true,data_initial.parameterization)

    para_MCMC.N = data_initial1.N_MCMC;
    start_point = parametrization(D_ini_only,data_initial.parameterization);
    result_emu_initial_only = MCMC_tensor_emulator(data_initial1,para_MCMC,start_point);
    result_emu_initial_only_h = MCMC_tensor_emulator_hierarchical(data_initial1,para_MCMC,start_point);

    RESULT.result_emu_initial=result_emu_initial;
    RESULT.result_emu_initial_h=result_emu_initial_h;
    RESULT.D_ini=D_ini;
    RESULT.hessian_ini=hessian_ini;
    RESULT.D_ini_h=D_ini_h;
    RESULT.hessian_ini_h=hessian_ini_h;

    RESULT.result_emu_seq=result_emu_seq;
    RESULT.D_seq=D_seq;
    RESULT.hessian_seq=hessian_seq;
    RESULT.data_seq=data_seq;

    RESULT.result_emu_seq_lik=result_emu_seq_lik;
    RESULT.D_seq_lik=D_seq_lik;
    RESULT.hessian_seq_lik=hessian_seq_lik;
    RESULT.data_seq_lik=data_seq_lik;

    RESULT.result_emu_initial_only=result_emu_initial_only;
    RESULT.result_emu_initial_only_h=result_emu_initial_only_h;
    RESULT.D_ini_only=D_ini_only;
    RESULT.D_ini_only_h=D_ini_only_h;
    RESULT.hessian_ini_only=hessian_ini_only;
    RESULT.hessian_ini_only_h=hessian_ini_only_h;

    RESULT.D_true=D_true;
    RESULT.Err=Err;
end