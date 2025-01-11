%% this file reads best parameters from optimization and plot vmax and energy 
clc
clear
close all

%% 
fn = 'purekinetic_Sept_2024_extrarun_highera.txt';
A = readmatrix(fn);
[M,I] = min(A(:,end));
pars = A(I,:);
pars = pars(1:end-1);

%% concentration  
TF_conc_t = get_aps(2:9);
TF_conc_const = zeros(1,8)+TF_conc_t(1);
nbd = 5;
nbd_s = 4;
nbd_v = 4;
mut_mat_s = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_v = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_sv = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];

load('new_promoter_activity_single.mat');
% load('random_dataset.mat');

real_data_Ps = new_real_data_Ps(1:end-1,:); % exclude the real WT, instead WT is Ps only, no mutations
% real_data_Ps = Ps_rand;
real_data_Ps_std = new_Ps_only_std(1:end-1,:);

real_data_Pv = new_real_data_Pv(1:end-1,:);
% real_data_Pv = Pv_rand;
real_data_Pv_std = new_Pv_only_std(1:end-1,:);

load('new_PsPv_promoter_activity.mat')
real_data_PsPv = PvPs_avg_new(1:end-1,:);
real_data_Psv_std = PvPs_std_new(1:end-1,:);

%% vmax grouping
group_array_s.g1 = 1; 
group_array_s.g2 = 2; 
group_array_s.g3 = 3; 
group_array_s.g4 = 4; 
group_array_s.g5 = 5; 
group_array_s.g6 = 6; 
group_array_s.g7 = 7; 
group_array_s.g8 = 8; 

group_array_v.g1 = 1; 
group_array_v.g2 = 2; 
group_array_v.g3 = 3; 
group_array_v.g4 = 4; 
group_array_v.g5 = 5; 
group_array_v.g6 = 6; 
group_array_v.g7 = 7; 
group_array_v.g8 = 8; 

group_array_sv.g1 = 1; 
group_array_sv.g2 = 2; 
group_array_sv.g3 = 3; 
group_array_sv.g4 = 4; 
group_array_sv.g5 = 5; 
group_array_sv.g6 = 6; 
group_array_sv.g7 = 7; 
group_array_sv.g8 = 8; 

% %% generate simulated data
% [final_diff,sim_data_s,sim_data_v,sim_data_sv,vmax_s_final,vmax_v_final,vmax_sv_final,vmax_array_s,vmax_array_v,vmax_array_sv,...
%     energyi_s,energyi_v,energyi_sv,H_conc_t,A_conc_t] = objective_function(...
%     nbd,pars,TF_conc_t,...
%     mut_mat_s,mut_mat_v,mut_mat_sv, ...
%     real_data_Ps,real_data_Pv,real_data_PsPv,...
%     group_array_s,group_array_v,group_array_sv);
% save('sim_data_original.mat','sim_data_s','sim_data_v','sim_data_sv')
% save('best_fit_energy.mat','energyi_s','energyi_sv','energyi_v')
% save('best_fit_vmax.mat','vmax_s_final','vmax_v_final','vmax_array_s','vmax_array_v','vmax_array_sv','group_array_s','group_array_v')
% save('best_fit_RNAP_t.mat','H_conc_t','A_conc_t')

%%
% sim data with constant RNA pol
[~,sim_data_s_constsig,sim_data_v_constsig,sim_data_sv_constsig,~,~,~,~,~,~,...
    ~,~,~,H_conc_t_norm,A_conc_t_norm] = objective_function_constant_pol(...
    nbd,pars,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv);
save('sim_data_const_holoenzyme.mat','sim_data_s_constsig','sim_data_v_constsig','sim_data_sv_constsig')

%%
Comparing_prediction_btwstrains_constant_pol