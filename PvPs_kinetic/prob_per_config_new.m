function p = prob_per_config_new(nbd, config, energyi,mut,TF_conc)

BF = Boltzmann_factor_per_config(config,energyi,TF_conc);
% disp(BF)
Z = partition_function(nbd,energyi,mut,TF_conc);
% disp(Z)

p = BF/Z;

end 

