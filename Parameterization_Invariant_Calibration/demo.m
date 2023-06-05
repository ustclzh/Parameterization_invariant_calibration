addpath(genpath(pwd));

distance_gp
inst
labels = {'uni+uni','matf+matf','wish(I,d+1)+wish','matf+uni','wish+uni'};
DATA = {data1,data2,data4,data12,data14};
%DATA = {data_initial11,data_initial14};
LABEL = labels;
k = length(DATA); 
Theta_true = DATA{3}.THETA_true;
for i = 1:60
    err = [];
    D_true=invparametrization(Theta_true(i,:))
       temp_time = tic;
    a = randn(d,rem(temp_time,1000));
    a = a*a';
    [s,~,~]=svd(a);
    D_true=s*diag(True_eig(1:d))*s';
    data_new_true=Model_Generator(D_true,data1); % data
    for j =1 : k
        data_current=Model_Generator(D_true,DATA{j},data_new_true); % data
        RESULT{j}=Posterior_Infrerence(data_current, D_true, strategy, 1, design_criterion);
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