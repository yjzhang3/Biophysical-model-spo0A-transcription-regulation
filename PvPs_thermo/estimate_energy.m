function [pars,M] = estimate_energy(...
    group_array_s,group_array_v,group_array_sv,file,ID,maxtime,strain_type)
% main script for solving parameters

% strain_type = [2,4,6];

load('new_promoter_activity_single.mat');
load('new_PsPv_promoter_activity.mat');

real_data_Ps = new_real_data_Ps(1:end-1,strain_type); % exclude the real WT, instead WT is Ps only, no mutations
% real_data_Ps = Ps_rand;
% real_data_Ps_std = new_Ps_only_std(1:end-1,:);

real_data_Pv = new_real_data_Pv(1:end-1,strain_type);
% real_data_Pv = Pv_rand;
% real_data_Pv_std = new_Pv_only_std(1:end-1,:);

real_data_PsPv = PvPs_avg_new(1:end-1,strain_type);


%% assign parameters
TF_conc_t = get_aps(2:9);
nbd = 5;
mut_mat_s = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_v = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_sv = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_s = mut_mat_s(strain_type,:);
mut_mat_v = mut_mat_v(strain_type,:);
mut_mat_sv = mut_mat_sv(strain_type,:);

%% set up bounds 

% 1,2,3,s,v, (1-5)
% 12,13,1s,1v, (6-9)
% 23,2s,2v, (10-12)
% 3s,3v, (13,14)
% sv (15)
% 123,12s,12v (16-18)
% 13s,13v (19-20)
% 1sv (21)
% 23s, 23v (22-23)
% 2sv (24)
% 3sv (25)
% 123s 123v (26-27)
% 12sv 13sv (28-29)
% 23sv (30)
ne_raw = nbd+nchoosek(nbd,2)+nchoosek(nbd,3)+nchoosek(nbd,4);
ne = ne_raw-length(mut_mat_s(1,:))+1;
nvars = ne+numel(fieldnames(group_array_s))+numel(fieldnames(group_array_v))+2*length(TF_conc_t); 
% nvars = ne;  

lb_overall = zeros(1,nvars)+exp(-5);
ub_overall = zeros(1,nvars)+exp(5);

lb_overall(1) = 0.2;
ub_overall(1) = 1; % alpha can't be greater than 1 (pars1 is alpha directly)


lb_overall(2:3) = 1;
ub_overall(2:3) = 1; % keep promoter energy fixed because we know it's embedded in concentration terms


lb_overall([15-2,21-2,24-2,25-2,28-2,29-2,30-2]) = 1;
ub_overall([15-2,21-2,24-2,25-2,28-2,29-2,30-2]) = 1; % no interaction between promoters

% % vmax 
fng_s = fieldnames(group_array_s);
fng_v = fieldnames(group_array_v);
ns = numel(fng_s);
nv = numel(fng_v);
ntot = ns+nv;
% 
lb_overall(ne+1:ne+ntot) = 400;
ub_overall(ne+1:ne+ntot) = 900;

% % RNAp concentration 
lb_overall(ne+ntot+1:end) = 0;
ub_overall(ne+ntot+1:end) = 2;

%% start to fit data
[pars,M] = fit_data(...
    nbd,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv,...
    nvars,lb_overall,ub_overall,file,ID,maxtime);

