function [diff,sim_data] = fit_Hill_diffK(p,real_data)

K_Pnotspo0A = p(1); % 1/exp(-G_promoter)
K_Pspo0A = p(2);

vmax_Pnotspo0A = p(3); % maximum initiation rate
R_time = p(4:10); % time-dependent RNAP

sim_Pnotspo0A = zeros(length(R_time),1);
sim_Pspo0A = zeros(length(R_time),1);
for i = 1:length(R_time)
    sim_Pnotspo0A(i) = vmax_Pnotspo0A*R_time(i)/(R_time(i) + K_Pnotspo0A);
    sim_Pspo0A(i) = vmax_Pnotspo0A*R_time(i)/(R_time(i) + K_Pspo0A);
end
sim_data = [sim_Pnotspo0A,sim_Pspo0A];

diff = 1/length(sim_data(:))*sum((real_data(:)-sim_data(:)).^2./real_data(:).^2);