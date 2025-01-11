%% table S7
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

