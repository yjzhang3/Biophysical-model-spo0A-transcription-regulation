function [onall,m] = prob_over_time_allboxbound_nopromoter(n_0Abox,energyi,mut,TF_conc_t)
% compute the probability of each configuration in a system with no
% promoter, and output the last row, which is probability of all boxes
% bound

bins = dec2bin(0:(2^n_0Abox-1), n_0Abox) - '0';

final_config = bins;
for ss = 1:length(mut) % ss range from 1 to 3 
    if mut(ss) == 0
        final_ind = find(final_config(:,ss)~=1); % final_ind stores the index 
        % with respect to 8 configurations
        final_config = final_config(final_ind,:);
    end
end

m = zeros(length(final_config(:,1)),length(TF_conc_t)); % each row is one configuration's time-dependent probability  (consistent with above that each configuration is each row)

for i = 1:length(final_config(:,1))
    curr_config = final_config(i,:);
    for j = 1:length(TF_conc_t)
        m(i,j) = prob_per_config_new(n_0Abox, curr_config, energyi,mut,TF_conc_t(j));% use the set of functions that do not depend on RNAP concentration
    end
end

onall = m(end,:); % probability that all  boxes are bound (as large as the nubmer of boxes)
