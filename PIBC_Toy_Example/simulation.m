addpath(genpath(pwd));
load("True_parameter2.mat");
load('design.mat')
d = 2;
D_true=D_TRUE{d,200};
True=pdeparameter(D_true); % data

for i =1:4
    DATA{i} = train_emulator(DATA2{i});
end 

% true parameter
for all_cycle = 1:100
D_true=D_TRUE{d,all_cycle+inst*30};
D_true
N=1;
True=pdeparameter(D_true); % data
for i =1:4
DATA{i}.obs=True.obs;
end
% optimization based on simulator
%RESULT_SIMU=run_simulator_based_infer(data_initial);
% optimization based on emulator+initial desin
RESULT_SP{all_cycle}=run_emulator_based_infer(DATA{1},DATA{2},D_true);
RESULT_FR{all_cycle}=run_emulator_based_infer(DATA{3},DATA{4},D_true);
RESULT_SIM{all_cycle}=run_simulator_based_infer(DATA{1});

filename = ['eg1mcmc','instance',num2str(inst),'.mat'];
save(filename)
end