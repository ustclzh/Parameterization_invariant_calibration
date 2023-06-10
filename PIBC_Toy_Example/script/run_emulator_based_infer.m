function RESULT=run_emulator_based_infer(data_initial,data_initial1,D_true)
    %%
        [D_ini,x_ini,hessian_ini]=map_optimization(data_initial);    
        para_MCMC.N = data_initial.N_MCMC;
        start_point = parametrization(D_ini,data_initial.parameterization);
        result_emu_initial = MCMC_tensor_emulator(data_initial,para_MCMC,start_point);
        disp('MAP of initial design based Posterior and err')
        D_ini_MCMC=invparametrization(MAP_est(result_emu_initial.chain(1000:end,:)));
        D_ini
        Err(1)=distance_matrix(D_ini,D_true,data_initial.parameterization)
        Err(4)=distance_matrix(D_ini_MCMC,D_true,data_initial.parameterization)
        % optimization based on emulator+initial desin+follow up design
        [D_seq,x_seq,hessian_seq,data_seq]=sequential_design(data_initial);
        para_MCMC.N = data_initial.N_MCMC;
        start_point = parametrization(D_seq,data_initial.parameterization);
        result_emu_seq = MCMC_tensor_emulator(data_seq,para_MCMC,start_point);
        disp('MAP of sequential design based Posterior and err')
        D_seq_MCMC=invparametrization(MAP_est(result_emu_seq.chain(1000:end,:)));
        D_seq 
        Err(2)=distance_matrix(D_seq,D_true,data_initial.parameterization)
        Err(5)=distance_matrix(D_seq_MCMC,D_true,data_initial.parameterization)
        % optimization based on emulator+support point+follow up design
        [D_ini_only,x_ini_only,hessian_ini_only]=map_optimization(data_initial1);
        para_MCMC.N = data_initial1.N_MCMC;
        start_point = parametrization(D_ini_only,data_initial.parameterization);
        result_emu_initial_only = MCMC_tensor_emulator(data_initial1,para_MCMC,start_point);
        D_ini_only_MCMC=invparametrization(MAP_est(result_emu_initial_only.chain(1000:end,:)));
        disp('MAP of SP design based Posterior and err')
        D_ini_only 
        Err(3)=distance_matrix(D_ini_only,D_true,data_initial.parameterization)
        Err(6)=distance_matrix(D_ini_only_MCMC,D_true,data_initial.parameterization)
    
        RESULT.result_emu_initial=result_emu_initial;
        RESULT.D_ini=D_ini;
        RESULT.D_ini_MCMC=D_ini_MCMC;
        RESULT.hessian_ini=hessian_ini;
    
    
        RESULT.result_emu_seq=result_emu_seq;
        RESULT.D_seq=D_seq;
        RESULT.D_seq_MCMC=D_seq_MCMC;
        RESULT.hessian_seq=hessian_seq;
        RESULT.data_seq=data_seq;
        
        RESULT.result_emu_initial_only=result_emu_initial_only;
        RESULT.D_ini_only=D_ini_only;
        RESULT.D_ini_only_MCMC=D_ini_only_MCMC;
        RESULT.hessian_ini_only=hessian_ini_only;
    
        RESULT.D_true=D_true;
        RESULT.Err=Err;
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
    
    
    
    