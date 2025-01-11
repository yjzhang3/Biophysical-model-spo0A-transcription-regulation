% examine the occupancy of 0A boxes without consider promoter or vmax
%% use best-fit energy to extract energy of 0A only (no promoter)
%%%%%%%%%Ps
% 1,2,3,s (1-4)
% 12,13,1s (5-7)
% 23,2s (8-9)
% 3s (10)
% 123, 12s (11-12)
% 13s (13)
% 23s (14)
% 123s (15)
clear
clc
close all
load('best_fit_energy.mat')
load('best_fit_model_prediction_colors.mat')
load('best_fit_vmax.mat')
energyi_nopromoter = energyi_s([1,2,3,5,6,8,11]);
TF_conc_t = get_aps(2:9);
% mut_mat_s = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_s_order_presence = [[1,0,0];[0,1,0];[0,0,1];[1,1,0];[1,0,1];[0,1,1];[1,1,1]]; % reorder in the order of presence, '1','2','3','4'...
%% compute the probability that all three 0A boxes are bound over time
clc
prob_all_bound_over_time_all_strain = zeros(length(TF_conc_t),length(mut_mat_s_order_presence(:,1))); % 8 rows for 8 time points, 8 columns for 8 strains, 

for si = 1:length(mut_mat_s_order_presence(:,1))  
    [onall,m] = prob_over_time_allboxbound_nopromoter(3,energyi_nopromoter,mut_mat_s_order_presence(si,:),TF_conc_t);
    prob_all_bound_over_time_all_strain(:,si) = onall;
end

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
mycolors = hex2rgb(hexColor); % order of color
plot_all_sim(prob_all_bound_over_time_all_strain,flipud(mycolors([8,1:7],:))) % sim data must have each strain data on each column
legend('1','2','3','12','13','23','123')
xlabel('Time (h)')
xlim([2 9])
ylabel('Probability')
set(gca,'FontSize',28)
% title('Probability that all 0A boxes are bound')
f = gcf;
exportgraphics(f,'best_fit_allboxbound_nopromoter.tif','Resolution',600)