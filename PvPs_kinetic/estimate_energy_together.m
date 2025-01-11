function [pars,M] = estimate_energy_together(...
    group_array_s,group_array_v,group_array_sv,filename,ID,maxtime)
% this finds both the energy values and the time-dependent [RNAPH] for
% promoter Ps because we would like to keep the number of paramters smaller
% than the number of real data points with 0A4


%% load data and choose strains
% *********** loading the right dataset? ***********
% original
% random
% extra

load('new_promoter_activity_single.mat');
load('new_PsPv_promoter_activity.mat')

real_data_Ps = new_real_data_Ps(1:end-1,:); % exclude the real WT, instead WT is Ps only, no mutations
% real_data_Ps = Ps_rand;
% real_data_Ps_std = new_Ps_only_std(1:end-1,:);

real_data_Pv = new_real_data_Pv(1:end-1,:);
% real_data_Pv = Pv_rand;
% real_data_Pv_std = new_Pv_only_std(1:end-1,:);

real_data_PsPv = PvPs_avg_new(1:end-1,:);
% real_data_Psv_std = PvPs_std_new(1:end-1,:);


%% assign known parameters
% *********** number of binding sites for each promoter correct? ***********
% 4 for Pv, 5 for Ps
% 4 for both Pv and Ps
TF_conc_t = get_aps(2:9);
nbd = 5;
mut_mat_s = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_v = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_sv = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];

%% set up bounds for unknowns

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
nvars = ne+8+8+2*length(TF_conc_t);

lb_overall = zeros(1,nvars)+exp(-3);
ub_overall = zeros(1,nvars)+exp(14);

lb_overall(1) = 0.2;
ub_overall(1) = 1; % alpha can't be greater than 1 (pars1 is alpha directly)

lb_overall(2:3) = 1;
ub_overall(2:3) = 1; % keep promoter energy fixed because we know it's embedded in concentration terms

% 1 2 3 12 13 23 123
% X X X  6 7 10 16
% load('aug_solved_single_site_energy.mat')
% % x1 = min(pars_range); %lb
% % x2 = max(pars_range); %ub
% lb_overall([6-2]) = 1;
% ub_overall([6-2]) = 1;
% lb_overall([16-2]) = exp(-x1(7));
% ub_overall([16-2]) = exp(-x2(7));
% lb_overall([6-2]) = exp(-x1(4));
% ub_overall([6-2]) = exp(-x2(4));
% lb_overall([7-2]) = exp(-x2(5));
% ub_overall([7-2]) = exp(-x1(5));

lb_overall([15-2,21-2,24-2,25-2,28-2,29-2,30-2]) = 1;
ub_overall([15-2,21-2,24-2,25-2,28-2,29-2,30-2]) = 1; % no interaction between promoters

lb_overall([8-2,11-2,13-2,17-2,19-2,22-2,26-2]) = 1;
ub_overall([8-2,11-2,13-2,17-2,19-2,22-2,26-2]) = 1; % no TF-polH interaction

lb_overall([9-2,12-2,14-2,18-2,20-2,23-2,27-2]) = 1;
ub_overall([9-2,12-2,14-2,18-2,20-2,23-2,27-2]) = 1; % no TF-polA interaction 


% vmax should be positive
lb_overall(ne+1:ne+16) = 10;
ub_overall(ne+1:ne+16) = 1500;

% RNApH and RNApA
lb_overall(ne+16+1:end) = 0;
ub_overall(ne+16+1:end) = 2;

%% fit to real data
% *********** is standard deviation part of the objective function? ***********
% Yes
% No
[pars,M] = fit_data_together_wpenalty(...
    nbd,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv,real_data_Ps,real_data_Pv,real_data_PsPv, ...
    lb_overall,ub_overall,...
    nvars,group_array_s,group_array_v,group_array_sv,filename,ID,maxtime);

