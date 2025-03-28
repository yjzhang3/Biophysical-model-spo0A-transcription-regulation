function [pars,M] = estimate_energy(...
    group_array_s,group_array_v,group_array_sv,file,ID)
% main script for solving parameters

load('new_promoter_activity_single.mat');
load('new_PsPv_promoter_activity.mat');

real_data_Ps = new_real_data_Ps(1:end-1,:); % exclude the real WT, instead WT is Ps only, no mutations
% real_data_Ps = Ps_rand;
% real_data_Ps_std = new_Ps_only_std(1:end-1,:);

real_data_Pv = new_real_data_Pv(1:end-1,:);
% real_data_Pv = Pv_rand;
% real_data_Pv_std = new_Pv_only_std(1:end-1,:);

real_data_PsPv = PvPs_avg_new(1:end-1,:);


%% assign parameters
TF_conc_t = get_aps(2:9);
nbd = 5;
mut_mat_s = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_v = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_sv = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];

%% set up bounds 

% energy related pars
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

ne = nbd+nchoosek(nbd,2)+nchoosek(nbd,3)+nchoosek(nbd,4)-3+1;
% nvars = ne+7+length(TF_conc_t); 
nvars = ne+8+numel(fieldnames(group_array_v))+2*length(TF_conc_t); 
% nvars = ne+7+numel(fieldnames(group_array_v))+length(TF_conc_t); 

lb_overall = zeros(1,nvars)-3;
ub_overall = zeros(1,nvars)+8;

lb_overall(1) = 0.2;
ub_overall(1) = 0.8;

lb_overall([2:3]) = 0;
ub_overall([2:3]) = 0;

lb_overall([15-2,21-2,24-2,25-2,28-2,29-2,30-2]) = 0;
ub_overall([15-2,21-2,24-2,25-2,28-2,29-2,30-2]) = 0; % no interaction between promoters

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

% Ps is kinetic (1s,2s,3s,12s,13s,23s,123s = 0, otherwise it's not kinetic)
lb_overall([8-2,11-2,13-2,17-2,19-2,22-2,26-2]) = 0;
ub_overall([8-2,11-2,13-2,17-2,19-2,22-2,26-2]) = 0;

% vmax 
fng_v = fieldnames(group_array_v);
nv = numel(fng_v); %1
ns = 8;
ntot = ns+nv;
% ntot = ns;

lb_overall(ne+1:ne+ntot) = 10;
ub_overall(ne+1:ne+ntot) = 2000;

% RNApH and RNApA concentration 
lb_overall(ne+ntot+1:ne+ntot+2*length(TF_conc_t)) = 0;
ub_overall(ne+ntot+1:ne+ntot+2*length(TF_conc_t)) = 2;
% lb_overall(ne+ntot+1:ne+ntot+length(TF_conc_t)) = 0;
% ub_overall(ne+ntot+1:ne+ntot+length(TF_conc_t)) = 2;


%% start to fit data
[pars,M] = fit_data(...
    nbd,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv,...
    nvars,lb_overall,ub_overall,file,ID);
end