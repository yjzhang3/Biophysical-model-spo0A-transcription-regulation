function vmax_final = vmax_assign(n_conf,vmax_array,group_array)
% this function helps us determine how many unknown vmax we want to have
% (for example), some configurations can have the same vmax

% input: number of configurations, a struct that specifies vmax for each group of configurations
% and a struct that specifies which strains should share the same vmax

% output: a vector of vmax with the same length as the total number of
% configurations (allow repetitions, because eseentially each position
% corresponds to each configuration)

vmax_final = zeros(n_conf,1);

fn = fieldnames(vmax_array);
for k=1:numel(fn)
    curr_val = vmax_array.(fn{k});
    curr_ind = group_array.(fn{k});
    vmax_final(curr_ind) = curr_val;
end
end
