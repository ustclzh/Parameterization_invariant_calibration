% note: true value generated from the prior model. 
% compare the misspecified prior and correctly specified prior
% compare design 
addpath(genpath(pwd));
%eg = 3;
%distance_gp = 5;
%ini_size = 5;
%seq_size = 5;
%inst = 1;
%design_algo = 2;
%True_theta = 1;
%design_algo = 2;
%design_criterion = 2
strategy = 1;
do_mcmc = 1;

data=config(eg,distance_gp,ini_size,seq_size);% eg, decomposition
d = data.d;
D_true = eye(d);
data.iso_kernel=1;
data.design_algo = design_algo;
data0=Model_Generator(D_true,data); % data

data1=Design_Initial(data0,1);
data12=Design_Initial(data0,1,2);
data2=Design_Initial(data0,2);
%data03=Design_Initial(data0,3);
data14=Design_Initial(data0,1,4, eye(d), d+1);
data4=Design_Initial(data0,4,4, eye(d), d+1);
data1=Emulator_Train(data1);
data2=Emulator_Train(data2);
data12=Emulator_Train(data12);
%data_iso_3=Emulator_Train(data03);
data4=Emulator_Train(data4);
data14=Emulator_Train(data14);

distance_gp
inst
labels = {'uni+uni','matf+matf','wish(I,d+1)+wish','matf+uni','wish+uni'};
DATA = {data1,data2,data4,data12,data14};
%DATA = {data_initial11,data_initial14};
LABEL = labels;
k = length(DATA); 
Theta_true = DATA{True_theta}.THETA_true;
for i = 1:60
    err = [];
    D_true=invparametrization(Theta_true(i,:));
    data_new_true=Model_Generator(D_true,data1); % data
    for j =1 : k
        data_current=Model_Generator(D_true,DATA{j},data_new_true); % data
        RESULT{j}=Posterior_Infrerence(data_current, D_true, strategy, do_mcmc, design_criterion);
%         if j ==1 || j == 3
%             RESULT{j}.bo2 = bayesian_optimization_2(data_current,3, 30);
%         end
        err(j,:) = [RESULT{j}.err_ini, RESULT{j}.err_seq];
    end
    %err = [err; RESULT{1}.bo2.RMSE(end), RESULT{3}.bo2.RMSE(end)]
    err
    %
    ERR(i,:) = err(:)';
    RES_prop{i} = RESULT;
    filename = ['res/eg',num2str(eg),'ini',num2str(ini_size),'seq',num2str(seq_size),...
        'design_alg',num2str(design_algo),'True_theta',num2str(True_theta),'GP_type',...
        num2str(distance_gp),'design_criterion',num2str(design_criterion),'instance',num2str(inst),'.mat'];
    save(filename)
end


