function BF = Boltzmann_factor_per_config(config,energyi,TF_conc)
% only calculate the bolzmann factor as if promoter is not there

% G0 = stand_energy(config,energyi); % we want to try optimizing directly e(-G) term, not G
G0 = stand_energy_linear(config,energyi);

num_TF = sum(config); % number of TF molecules (configurations are initialized without promoter !)

conc_TF = TF_conc;

BF = G0*(conc_TF)^(4*num_TF);

end