function [RMSE,y_pred,y_simu,D_test] = Compare_GP(data0,D_test)
if nargin ==1
    data = data0;
    design_para = data.design_para;
    design_matrix = data.design_matrix;
    n = size(design_para,1);
    n_predict = round(n*0.2);
    n_train = n-n_predict;
    ind = randsample(n,n);
    ind_test = ind(1:n_predict);
    ind_train = ind(n_predict+1:n);
    data.design_para = design_para(ind_train,:);
    data.design_matrix = data.design_matrix(:,:,ind_train);
    y_simu = data.Y_simulator(:,ind_test)';
    data.Y_simulator = data.Y_simulator(:,ind_train);
    data.q_0 = size(data.design_para,1);
    D_test = design_para(ind_test,:);
    n_test = size(D_test,1);
else
    n_test = size(D_test,1);
    y_simu = zeros(n_test,length(data.obs));
    for i =1 : n_test
        y_simu(i,:) = simulator(data,D_test(i,:));
    end

end
data = Emulator_Train(data,1);
y_pred = zeros(n_test,length(data.obs));
for i = 1:n_test
y_pred(i,:) = Emulator_Predict(data,D_test(i,:));
end
RMSE = (mean(mean((y_simu-y_pred).^2)).^(1/2));
end