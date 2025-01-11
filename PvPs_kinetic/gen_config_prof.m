function [n_arr,c_arr] = gen_config_prof(config,TF_conc,RNAp_conc)
% input: a specific configuration (specified by 01 array)
% output: the number and concentration profile for each configuration

n_arr = zeros(2,1); % assume only two types of protein, TF and RNAP
c_arr = zeros(2,1); % assume only two types of protein, TF and RNAP

if config(end) == 1 % assume RNAP is of one molecule
    n_arr(end) = 1;
end

c_arr(end) = RNAp_conc; 

% assume besides RNAP, all the other binding sites are bound by *one* type
% of TF
n_arr(1) = sum(config(1:end-1)); 
c_arr(1) = TF_conc;

end
    