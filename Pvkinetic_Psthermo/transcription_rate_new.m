function TR = transcription_rate_new(nbd,energyi,mut,TF_conc,RNAp_conc,group_array,vmax_array)
% this is another way to incorporate mutations in the data

%% exclude configurations that have no promoters bound
bins = dec2bin(0:(2^nbd-1), nbd) - '0';
on_ind = find(bins(:,end) == 1); 
on_config = bins(on_ind,:);

% consider mutated sites (exclude those configurations where mutated site
% is 1)
final_config = on_config;
for ss = 1:length(mut) % ss range from 1 to 3 or 4
    if mut(ss) == 0
        final_ind = find(final_config(:,ss)~=1); % final_ind stores the index 
        % with respect to 8 configurations
        final_config = final_config(final_ind,:);
    end
end

final_ind = find(ismember(on_config,final_config,'row'));

n_conf = 2^(nbd-1); % number of possible configurations excludig the promoter state

% vmax is now an array
% create a map that maps a unique vmax to each configuration

% lets make vmax assignment using a function instead
vmax_final = vmax_assign(n_conf,vmax_array,group_array);
% the numbering of each configuration should be as if there are only 0A
% boxes present 


TR = 0;
for oo = 1:length(final_config(:,1))
    if isvector(final_config)
        config = final_config;
    else
        config = final_config(oo,:);
    end
    curr_vmax = vmax_final(final_ind(oo));
    curr_p = prob_per_config_new(nbd,config,energyi,mut,TF_conc,RNAp_conc);
    curr_tr = curr_vmax*curr_p;
    TR = TR + curr_tr;
end

end