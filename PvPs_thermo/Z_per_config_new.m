function Z = Z_per_config_new(config, energyi,TF_conc,RNAp_conc)

nbd = length(config);
if length(energyi) == nbd+nchoosek(nbd,2)
    G0 = stand_energy_linear(config,energyi);
elseif length(energyi) == nbd+nchoosek(nbd,2)+nchoosek(nbd,3)
    G0 = stand_energy_linear_new(config,energyi);
elseif length(energyi) == nbd+nchoosek(nbd,2)+nchoosek(nbd,3)+nchoosek(nbd,4)
    G0 = stand_energy_linear_new_q(config,energyi);
end
% G0 = stand_energy_linear_new(config,energyi);

[num_arr,conc_arr] = gen_config_prof(config,TF_conc,RNAp_conc);

num_TF = num_arr(1)*4;
% num_TF = num_arr(1);
num_RNAp = num_arr(end);

conc_TF = conc_arr(1);
conc_RNAp = conc_arr(end);


Z = G0*(conc_TF^num_TF)*(conc_RNAp^num_RNAp);

end