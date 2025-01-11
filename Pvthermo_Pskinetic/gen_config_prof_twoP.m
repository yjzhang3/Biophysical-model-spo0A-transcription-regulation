function [n_arr,c_arr] = gen_config_prof_twoP(config,TF_conc,RNApH_conc,RNApA_conc)
% input: a specific configuration (specified by 01 array)
% output: the number and concentration profile for each configuration

n_arr = zeros(3,1); % assume only 3 types of protein, TF, RNAP-sigmaH, RNAP-sigmaA
c_arr = zeros(3,1); % assume only 3 types of protein, TF, RNAP-sigmaH, RNAP-sigmaA

% Let's always put Ps before Pv
n_arr(end) = config(end); % Pv-sigmaA
c_arr(end) = RNApA_conc; 

n_arr(end-1)= config(end-1); % Ps-sigmaH
c_arr(end-1) = RNApH_conc;

% assume besides RNAP, all the other binding sites are bound by *one* type
% of TF
n_arr(1) = sum(config(1:end-2)); % TF
c_arr(1) = TF_conc;

end
    