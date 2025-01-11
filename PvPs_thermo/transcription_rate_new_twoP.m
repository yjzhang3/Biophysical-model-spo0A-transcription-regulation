function TR_all = ...
    transcription_rate_new_twoP( ...
    nbd,energyi,mut,TF_conc,...
    RNApH_conc,RNApA_conc,...
    vmax_array_s,group_array_s,...
    vmax_array_v,group_array_v,...
    vmax_array_sv,group_array_sv)
% calculate the total probability of Ps and Pv promoter transcribed, if
% they transcribe independently

%% exclude configurations that have neither promoters bound

bins = dec2bin(0:(2^nbd-1), nbd) - '0';

Ps_final = find(bins(:,end-1) == 1 & bins(:,end) == 0); % XXXX10
Pv_final = find(bins(:,end-1) == 0 & bins(:,end) == 1);% XXXX01
PsPv_final = find(bins(:,end-1) == 1 & bins(:,end) == 1); % XXXX11

Ps_on_config = bins(Ps_final,:);
Pv_on_config = bins(Pv_final,:);
PsPv_on_config = bins(PsPv_final,:);

%% consider mutated sites (exclude those configurations where mutated site
% is 1)
Ps_final_config = Ps_on_config;
for ss = 1:length(mut) % ss range from 1 to 3
    if mut(ss) == 0
        final_ind = find(Ps_final_config(:,ss)~=1); % final_ind stores the index 
        % with respect to 8 configurations
        Ps_final_config = Ps_final_config(final_ind,:);
    end
end
Ps_final_ind = find(ismember(Ps_on_config,Ps_final_config,'row'));

clear final_ind
Pv_final_config = Pv_on_config;
for ss = 1:length(mut) % ss range from 1 to 3
    if mut(ss) == 0
        final_ind = find(Pv_final_config(:,ss)~=1); % final_ind stores the index 
        % with respect to 8 configurations
        Pv_final_config = Pv_final_config(final_ind,:);
    end
end
Pv_final_ind = find(ismember(Pv_on_config,Pv_final_config,'row'));


clear final_ind
PsPv_final_config = PsPv_on_config;
for ss = 1:length(mut) % ss range from 1 to 3
    if mut(ss) == 0
        final_ind = find(PsPv_final_config(:,ss)~=1); % final_ind stores the index 
        % with respect to 8 configurations
        PsPv_final_config = PsPv_final_config(final_ind,:);
    end
end
PsPv_final_ind = find(ismember(PsPv_on_config,PsPv_final_config,'row'));


%% assign vmax to each configuration
n_conf = 2^(nbd-2); % number of possible configurations if excluding the state of promoter

vmax_final_Ps = vmax_assign(n_conf,vmax_array_s,group_array_s);
vmax_final_Pv = vmax_assign(n_conf,vmax_array_v,group_array_v);
vmax_final_PsPv = vmax_assign(n_conf,vmax_array_sv,group_array_sv);

%% calculate transcription rate for XXXX01, XXXX10, and XXXX11, respectively  
Pv_TR = 0;
Ps_TR = 0;
both_TR = 0;

% 1,2,3,4,s,v, (1-6)
% 12,13,14,1s,1v, (7-11)
% 23,24,2s,2v, (12-15)
% 34,3s,3v, (16-18)
% 4s,4v,sv (19-21)


for oo = 1:length(Ps_final_config(:,1))

    config = Ps_final_config(oo,:);

    vmax = vmax_final_Ps(Ps_final_ind(oo)); % weight that config with its vmax

    Ps_TR = Ps_TR + vmax*prob_per_config_new_twoP(nbd,config,energyi,mut,TF_conc,RNApH_conc,RNApA_conc);
end

clear oo
for oo = 1:length(Pv_final_config(:,1))

    config = Pv_final_config(oo,:);

    vmax = vmax_final_Pv(Pv_final_ind(oo)); % weight that config with its vmax

    Pv_TR = Pv_TR + vmax*prob_per_config_new_twoP(nbd,config,energyi,mut,TF_conc,RNApH_conc,RNApA_conc);
end

clear oo

for oo = 1:length(PsPv_final_config(:,1))

    config = PsPv_final_config(oo,:);

    vmax = vmax_final_PsPv(PsPv_final_ind(oo)); % weight that config with its vmax

    both_TR = both_TR + vmax*prob_per_config_new_twoP(nbd,config,energyi,mut,TF_conc,RNApH_conc,RNApA_conc);
end

TR_all = Ps_TR+Pv_TR+both_TR;

end