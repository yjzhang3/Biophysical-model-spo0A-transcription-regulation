%% objective function for estimating the single binding site energy 
function [diff,sim_data,sim_plot_data] = objective_G_together_strain_weight(n,G,TF_conc,real_data,n_exp_per_strain)
% estimate the single binding energy (and possibly interaction energy)
% altogether
% real data = fraction bound. Each column should be the spo0A dependent
% fraction bound for a particular strain. N columns means N strains
% if there are n expriemnts per strain, then solve for energy n times for
% that strain

G1 = G(1);
G2 = G(2);
G3 = G(3);

% sim_data = zeros(size(real_data));

sim_data_k1 = zeros(length(TF_conc),n_exp_per_strain(1)); % generate model data points that are as many as the real data points (each column is a 
% concentration dependent probabilty
for k = 1:n_exp_per_strain(1)
    sim_data_i = zeros(length(TF_conc),1);
    for i = 1:length(TF_conc)
        sim_data_i(i) = 1-1/(1+TF_conc(i)^n*exp(-G1)); % 0A1
    end
    sim_data_k1(:,k) = sim_data_i;
end

sim_data_k2 = zeros(length(TF_conc),n_exp_per_strain(2));
for k = 1:n_exp_per_strain(2)
    sim_data_i = zeros(length(TF_conc),1);
    for i = 1:length(TF_conc)
        sim_data_i(i) = 1-1/(1+TF_conc(i)^n*exp(-G2)); % 0A2
    end
    sim_data_k2(:,k) = sim_data_i;
end 

sim_data_k3 = zeros(length(TF_conc),n_exp_per_strain(3));
for k = 1:n_exp_per_strain(3)
    sim_data_i = zeros(length(TF_conc),1);
    for i = 1:length(TF_conc)
        sim_data_i(i) = 1-1/(1+TF_conc(i)^n*exp(-G3)); % 0A3
    end
    sim_data_k3(:,k) = sim_data_i;
end 

sim_data = [sim_data_k1,sim_data_k2,sim_data_k3];

% diff = weighted_msd(sim_data(:),real_data(:));

diff = 1/length(sim_data(:))*sum((sim_data(:)-real_data(:)).^2);

sim_plot_data = [mean(sim_data_k1,2),mean(sim_data_k2,2),mean(sim_data_k3,2)]; % for plotting purposes, only
% plot the average across n data sets (which is the same under the same G
% come on...)


end