%% increase in gene expression with constant 0A~P or constant RNA pol
clear
best_fit_model_prediction_constant_0AP

%%
close all
clear
best_fit_model_prediction_constant_pol
%%
clear
clc
close all
load('sim_data_const_holoenzyme.mat') % pre-run and saved
load('sim_data_const_0AP.mat')
%% constant 0A~P (increase due to sigma only)
fc_const_TF_s = sim_data_s_constTF./sim_data_s_constTF(1,:); 
fc_const_TF_v = sim_data_v_constTF./sim_data_v_constTF(1,:);
fc_const_TF_sv = sim_data_sv_constTF./sim_data_sv_constTF(1,:);

%% constant holoenzyme (increase due to 0A~P only)
fc_const_holoenzyme_s = sim_data_s_constsig./sim_data_s_constsig(1,:); 
fc_const_holoenzyme_v = sim_data_v_constsig./sim_data_v_constsig(1,:);
fc_const_holoenzyme_sv = sim_data_sv_constsig./sim_data_sv_constsig(1,:);

%% percent change in gene expression under constant holoenzyme 
pc_const_holoenzyme_s = (sim_data_s_constsig(end,:)-sim_data_s_constsig(end,end))./sim_data_s_constsig(end,end);
pc_const_holoenzyme_v = (sim_data_v_constsig(end,:)-sim_data_v_constsig(end,end))./sim_data_v_constsig(end,end);
pc_const_holoenzyme_sv = (sim_data_sv_constsig(end,:)-sim_data_sv_constsig(end,end))./sim_data_sv_constsig(end,end);
