function [final_diff,sim_data_Ps,sim_data_Pv,sim_data_PsPv,...
    vmax_s_final,vmax_v_final,vmax_sv_final,H_conc_t,A_conc_t,...
    energyi_s,energyi_v,energyi_sv] = objective_function(...
    nbd,p,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv)

%% not account for 0A4 but there are triple interactions
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

% without constraint
% ne = ne_raw;
% energyi_sv = p(1:ne);

% with constraint
energyi_sv = zeros(1,ne_raw);
% load('fixed_energy_interval_n4_111.mat')
% bindingE = mean([pars_kd_lb(1:3);pars_kd_ub(1:3)],1);
load('aug_energy_Dec2024.mat')
bindingE = mean(pars_range(:,1:3));

ne = ne_raw-length(bindingE)+1; % number of unknown pars related to energy
alpha = p(1);
beta = 4*log(alpha);
energyi_sv(1:3) = exp(-bindingE)*exp(-beta);
energyi_sv(4:end) = p(1+1:ne); 

%%%%%%% ps
% 1,2,3,s (1-4)
% 12,13,1s (5-7)
% 23,2s (8-9)
% 3s (10)
% 123, 12s (11-12)
% 13s (13)
% 23s (14)
% without constraint
% energyi_s = zeros(1,nbd-1+nchoosek(nbd-1,2)+nchoosek(nbd-1,3));
% energyi_s(1:4) = p([1:4]);
% energyi_s(5:7) = p([6:8]);
% energyi_s(8:9) = p([10,11]);
% energyi_s(10) = p(13);
% energyi_s(11:14) = p([16,17,19,22]);

% with constraint
energyi_s = zeros(1,nbd-1+nchoosek(nbd-1,2)+nchoosek(nbd-1,3)+nchoosek(nbd-1,4));
energyi_s(1:4) = energyi_sv([1:4]);
energyi_s(5:7) = energyi_sv([6:8]);
energyi_s(8:9) = energyi_sv([10,11]);
energyi_s(10) = energyi_sv(13);
energyi_s(11:14) = energyi_sv([16,17,19,22]);
energyi_s(15) = energyi_sv(26);

%%%%%%%% pv
% 1,2,3,v (1-4)
% 12,13,1v (5-7)
% 23,2v (8-9)
% 3v (10)
% 123, 12v (11-12)
% 13v (13)
% 23v (14)
% without constraint
% energyi_v = zeros(1,nbd-1+nchoosek(nbd-1,2)+nchoosek(nbd-1,3));
% energyi_v(1:4) = p([1:3,5]);
% energyi_v(5:7) = p([6,7,9]);
% energyi_v(8:9) = p([10,12]);
% energyi_v(10) = p(14);
% energyi_v(11:14) = p([16,18,20,23]);

% with constraint
energyi_v = zeros(1,nbd-1+nchoosek(nbd-1,2)+nchoosek(nbd-1,3)+nchoosek(nbd-1,4));
energyi_v(1:4) = energyi_sv([1:3,5]);
energyi_v(5:7) = energyi_sv([6,7,9]);
energyi_v(8:9) = energyi_sv([10,12]);
energyi_v(10) = energyi_sv(14);
energyi_v(11:14) = energyi_sv([16,18,20,23]);
energyi_v(15) = energyi_sv(27);

%% assign vmax (custom grouping)
% fng_s = fieldnames(group_array_s);
% fng_v = fieldnames(group_array_v);
% fng_sv = fieldnames(group_array_sv);
% ns = numel(fng_s);
% nv = numel(fng_v);
% % nsv = numel(fng_sv);
% ntot = ns+nv;
% 
% vmax_raw = p(ne+1:ne+ntot);
% vmax_s = vmax_raw(1:ns);
% vmax_v = vmax_raw((1+ns):(ns+nv));
% vmax_sv = zeros(1,length(vmax_v));
% for i = 1:length(vmax_v)
%     vmax_sv(i) = vmax_s(1)+vmax_v(i);
% end
% 
% for k=1:numel(fng_s) % iterate through k groups of group_array
%     vmax_array_s.(fng_s{k}) = vmax_s(k); % assign unknown parameters to that same name in vmax_array
% end
% 
% clear k
% for k=1:numel(fng_v) % iterate through k groups of group_array
%     vmax_array_v.(fng_v{k}) = vmax_v(k); % assign unknown parameters to that same name in vmax_array
% end
% 
% clear k
% for k=1:numel(fng_sv) % iterate through k groups of group_array
%     vmax_array_sv.(fng_sv{k}) = vmax_sv(k); % assign unknown parameters to that same name in vmax_array
% end

%% non-multiplicative vmax for Pv
% pv: v1,v2,v3,v12,v23,v123 (kinetic)
% ps: 1 vmax (thermo)

% load("vmax_rnap_Hill_constraint.mat",'Pv_Hill_Vmax_Gp')
fng_s = fieldnames(group_array_s);
fng_v = fieldnames(group_array_v);
fng_sv = fieldnames(group_array_sv);

%%%%% plug in Pv
ns = numel(fng_s);
nv = 8;
ntot = ns+nv;
vmax_raw = p(ne+1:ne+ntot);

vmax_s_final = vmax_raw(1:ns);
vmax_v_raw = vmax_raw(ns+1:ntot); % stores the optimized parameter 

% Pv: 000,001,010,011,100,101,110,111
% v1,v2,v3,v12,v23,v123 
vmax_v_final = [vmax_v_raw(8),vmax_v_raw(3),vmax_v_raw(2),vmax_v_raw(5),...
    vmax_v_raw(1),vmax_v_raw(7),vmax_v_raw(4),vmax_v_raw(6)];

%%%%%% plug in Ps
% ns = numel(fng_s);
% nv = 7;
% ntot = nv;
% % ntot = nv;
% vmax_raw = p(ne+1:ne+ntot);
% 
% vmax_s_final = Ps_Hill_Vmax_Gp(1);
% vmax_v_raw = vmax_raw(1:ntot); % stores the optimized parameter 
% 
% % Pv: 000,001,010,011,100,101,110,111
% % v1,v2,v3,v12,v23,v123 
% vmax_v_final = [vmax_v_raw(7),vmax_v_raw(3),vmax_v_raw(2),vmax_v_raw(5),...
%     vmax_v_raw(1),vmax_v_raw(1)*vmax_v_raw(3),vmax_v_raw(4),vmax_v_raw(6)];


vmax_sv_final = zeros(1,length(fng_v));

% every configuration has additive vamx from both ps and pv
for i = 1:length(fng_v)
    vmax_sv_final(i) = vmax_s_final(1)+vmax_v_final(i);
end

%% assign vmax to each configuration
for k=1:numel(fng_s) % iterate through k groups of group_array
    vmax_array_s.(fng_s{k}) = vmax_s_final(k); % assign unknown parameters to that same name in vmax_array
end
clear k
for k=1:numel(fng_v) % iterate through k groups of group_array
    vmax_array_v.(fng_v{k}) = vmax_v_final(k); % assign unknown parameters to that same name in vmax_array
end
clear k
for k=1:numel(fng_sv) % iterate through k groups of group_array
    vmax_array_sv.(fng_sv{k}) = vmax_sv_final(k); % assign unknown parameters to that same name in vmax_array
end

%% assign time-dependent sigmaH concentration

H_conc_t = p(ne+ntot+1:ne+ntot+length(TF_conc_t));
% H_conc_t = Ps_Hill_Vmax_Gp(2:9);
RNApH_conc_t = H_conc_t;

A_conc_t = p(ne+ntot+length(TF_conc_t)+1:ne+ntot+2*length(TF_conc_t));
% A_conc_t = p(ne+ntot+1:ne+ntot+length(TF_conc_t));
% A_conc_t = Pv_Hill_Vmax_Gp(2:9);
RNApA_conc_t = A_conc_t;

%% generate simulated data for all WT and mutated type

% all data follow this format:
% each row: every hour
% each column: every mutant

n = length(mut_mat_s(:,1));

%%%%%%%%%% Ps only
sim_data_Ps = zeros(length(TF_conc_t),n);

for mm = 1:n
    sim_data_Ps(:,mm) = time_dep_TR_new_wSigma( ...
        nbd-1,energyi_s,mut_mat_s(mm,:),TF_conc_t, ...
        RNApH_conc_t,group_array_s,vmax_array_s);
end


diff_s = weighted_msd(real_data_Ps(:),sim_data_Ps(:));

%%%%%%%%%% Pv only
sim_data_Pv = zeros(length(TF_conc_t),n);

clear mm
for mm = 1:n
    sim_data_Pv(:,mm) = time_dep_TR_new_wSigma( ...
        nbd-1,energyi_v,mut_mat_v(mm,:),TF_conc_t, ...
        RNApA_conc_t,group_array_v,vmax_array_v);
end


diff_v = weighted_msd(real_data_Pv(:),sim_data_Pv(:));


%%%%%%%%% PsPv together
sim_data_PsPv = zeros(length(TF_conc_t),n);

clear mm
for mm = 1:n
    sim_data_PsPv(:,mm) = time_dep_TR_new_wSigma_twoP( ...
    nbd,energyi_sv,mut_mat_sv(mm,:),TF_conc_t,...
    RNApH_conc_t,RNApA_conc_t,...
    vmax_array_s,group_array_s,...
    vmax_array_v,group_array_v,...
    vmax_array_sv,group_array_sv);
end

diff_sv = weighted_msd(real_data_PsPv(:),sim_data_PsPv(:));


%% Penalty
% 0A3 ~ 0A1 > 0A4 > 0A2
% penalty(ksmall,klarge,g1,g2)
% g1 = 0.1;
% g2 = 0.1;
% p1 = penalty(p(1),p(2),g1,g2); % 0A1 < 0A2
% p2 = penalty(p(3),p(2),g1,g2); % 0A3 < 0A2
% p4 = (p(3)-p(1))^2; % 0A3 ~ 0A1

final_diff = diff_s+diff_v+diff_sv;
% +p1+p2;
% +p4;

end