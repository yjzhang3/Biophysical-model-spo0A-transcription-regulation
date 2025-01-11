function Z = Z_per_config_new_twoP(config, energyi,TF_conc,RNApH_conc,RNApA_conc)

nbd = length(config);
if length(energyi) == nbd+nchoosek(nbd,2)
    G0 = stand_energy_linear(config,energyi);
elseif length(energyi) == nbd+nchoosek(nbd,2)+nchoosek(nbd,3)
    G0 = stand_energy_linear_new(config,energyi);
elseif length(energyi) == nbd+nchoosek(nbd,2)+nchoosek(nbd,3)+nchoosek(nbd,4)
    G0 = stand_energy_linear_new_q(config,energyi);
end
% G0 = stand_energy_linear_new(config,energyi);

[num_arr,conc_arr] = gen_config_prof_twoP(config,TF_conc,RNApH_conc,RNApA_conc);

num_TF = num_arr(1)*4;
% num_TF = num_arr(1);
num_RNApH = num_arr(2);
num_RNApA = num_arr(3);

conc_TF = conc_arr(1);
conc_RNApH = conc_arr(2);
conc_RNApA = conc_arr(3);

Z = G0*(conc_TF^num_TF)*(conc_RNApH^num_RNApH)*(conc_RNApA^num_RNApA);

end