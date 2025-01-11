function sim_data_multi_site = predict_multi0A_using_single_siteE(pars,TF_conc,n)
% predict binding probability of constructs with more than one 0A box,
% assuming no interaction energy.
G1 = pars(1);
G2 = pars(2);
G3 = pars(3);

sim_data_multi_site = zeros(length(TF_conc),4); % 3 columns for 0A12, 0A13, and 0A23, and 0A123; 11 rows for 11 TF concentrations

for i = 1:length(TF_conc)
    sim_data_multi_site(i,1) = 1-1/(1+TF_conc(i)^n*exp(-G1)+TF_conc(i)^n*exp(-G2)+TF_conc(i)^(2*n)*exp(-G1-G2)); % 0A12
end

for i = 1:length(TF_conc)
    sim_data_multi_site(i,2) = 1-1/(1+TF_conc(i)^n*exp(-G1)+TF_conc(i)^n*exp(-G3)+TF_conc(i)^(2*n)*exp(-G1-G3)); % 0A13
end

for i = 1:length(TF_conc)
    sim_data_multi_site(i,3) = 1-1/(1+TF_conc(i)^n*exp(-G2)+TF_conc(i)^n*exp(-G3)+TF_conc(i)^(2*n)*exp(-G2-G3)); % 0A23
end

for i = 1:length(TF_conc)
    sim_data_multi_site(i,4) = 1-1/(1+TF_conc(i)^n*exp(-G1)+TF_conc(i)^n*exp(-G2)+TF_conc(i)^n*exp(-G3)+...
        TF_conc(i)^(2*n)*exp(-G1-G2)+TF_conc(i)^(2*n)*exp(-G1-G3)+TF_conc(i)^(2*n)*exp(-G2-G3)+...
        TF_conc(i)^(3*n)*exp(-G1-G2-G3)); % 0A123
end
end
