clear
clc
close all
load('sigA_sigH_real_data.mat')
norm_sigA = Physpac_mean./Physpac_mean(1);
norm_sigH = citG_mean./citG_mean(1);
norm_sigA_std = Physpac_st./Physpac_mean(1);
norm_sigH_std = citG_st./citG_mean(1);

load('multiple_run_mean_job61000.mat')
norm_H_conc_t = H_conc_t_multiple_run_mean./H_conc_t_multiple_run_mean(1);
norm_H_conc_t_std = H_conc_t_multiple_run_std./H_conc_t_multiple_run_mean(1);
norm_A_conc_t = A_conc_t_multiple_run_mean./A_conc_t_multiple_run_mean(1);
norm_A_conc_t_std = A_conc_t_multiple_run_std./A_conc_t_multiple_run_mean(1);

%% plot sigma dynamics measured using promoters other than spo0A (no normalization)
clc
% colors = ['#0063dc';'#ff0084'];
colors = ['#FF0000';'#0063dc'];
mycolors = hex2rgb(colors);

figure('position',[188,100,1300/2,550])
t = 2:8;
plot_areaerrorbar_wosim(t,Physpac_mean,Physpac_st,mycolors(1,:))% error bar for real data
hold on
% plot_all_sim_2(t,norm_A_conc_t(1:end-1)',norm_A_conc_t_std(1:end-1)',mycolors(1,:)) % color patch for predicted data
% hold on
plot_areaerrorbar_wosim(t,citG_mean,citG_st,mycolors(2,:))% error bar for real data
% hold on
% plot_all_sim_2(t,norm_H_conc_t(1:end-1)',norm_H_conc_t_std(1:end-1)',mycolors(2,:)) % color patch for predicted data
% legend('Experimental \sigma^{A}','','Predicted \sigma^{A}','Experimental \sigma^{H}','','Predicted \sigma^{H}','location','northwest')
legend('Phy-spank (\sigma^{A})','PcitG (\sigma^{H})','location','northwest')

xlabel('Time (h)')
ylabel('Promoter Activity')
set(gca,'FontSize',22)

% xlim([2 9])
ylim([0 5600])
% legend('WT')
f = gcf;
exportgraphics(f,'SigmaDynamics_citGhyspank.tif','Resolution',600)

%% plot sigma dynamics measured using promoters other than spo0A (normalized by t = 2h)
clc
% colors = ['#0063dc';'#ff0084'];
colors = ['#FF0000';'#0063dc'];
mycolors = hex2rgb(colors);

figure('position',[188,100,1300/2,550])
t = 2:8;
plot_areaerrorbar_wosim(t,norm_sigA,norm_sigA_std,mycolors(1,:))% error bar for real data
hold on
% plot_all_sim_2(t,norm_A_conc_t(1:end-1)',norm_A_conc_t_std(1:end-1)',mycolors(1,:)) % color patch for predicted data
% hold on
plot_areaerrorbar_wosim(t,norm_sigH,norm_sigH_std,mycolors(2,:))% error bar for real data
% hold on
% plot_all_sim_2(t,norm_H_conc_t(1:end-1)',norm_H_conc_t_std(1:end-1)',mycolors(2,:)) % color patch for predicted data
% legend('Experimental \sigma^{A}','','Predicted \sigma^{A}','Experimental \sigma^{H}','','Predicted \sigma^{H}','location','northwest')
legend('Phy-spank (\sigma^{A})','PcitG (\sigma^{H})','location','northwest','FontSize',22)

xlabel('Time (h)')
ylabel('Fold Change')
set(gca,'FontSize',28)

% xlim([2 9])
% ylim([0 2.5])
% legend('WT')
f = gcf;
exportgraphics(f,'SigmaDynamics_citGhyspank_norm.tif','Resolution',600)

%% plot hy-spac and Pv normalized by t = 2h
colors = ['#FF0000';'#ffA000'];
mycolors = hex2rgb(colors);
figure('position',[188,100,1300/2,550])
plot_areaerrorbar_wosim(t,norm_sigA,norm_sigA_std,mycolors(1,:))% error bar for real data
hold on
plot_all_sim_2(t,norm_A_conc_t(1:end-1)',norm_A_conc_t_std(1:end-1)',mycolors(2,:)) % color patch for predicted data
hold on
xlabel('Time (h)')
ylabel('Fold Change')
set(gca,'FontSize',22)
legend('Phy-spank','','spo0A Pv','location','northwest')
f = gcf;
exportgraphics(f,'SigmaDynamics_hyspankspo0A_norm.tif','Resolution',600)

%% plot citG and Ps normalized by t = 2h
colors = ['#0063dc';'#00BfEB'];
mycolors = hex2rgb(colors);
figure('position',[188,100,1300/2,550])
plot_areaerrorbar_wosim(t,norm_sigH,norm_sigH_std,mycolors(1,:))% error bar for real data
hold on
plot_all_sim_2(t,norm_H_conc_t(1:end-1)',norm_H_conc_t_std(1:end-1)',mycolors(2,:)) % color patch for predicted data
hold on
xlabel('Time (h)')
ylabel('Fold Change')
set(gca,'FontSize',22)
legend('PcitG','','spo0A Ps','location','northwest')
f = gcf;
exportgraphics(f,'SigmaDynamics_citGspo0A_norm.tif','Resolution',600)

%% plot sigma dynamics measured using promoters of both spo0A and nonspo0A(no normalization)
clc
% colors = ['#0063dc';'#ff0084'];
colors = ['#FF0000';'#ffA000';'#0063dc';'#00BfEB'];
mycolors = hex2rgb(colors);

figure('position',[188,100,1300/2,550])
t = 2:8;
plot_areaerrorbar_wosim(t,norm_sigA,norm_sigA_std,mycolors(1,:))% error bar for real data
hold on
plot_all_sim_2(t,norm_A_conc_t(1:end-1)',norm_A_conc_t_std(1:end-1)',mycolors(2,:)) % color patch for predicted data
hold on
plot_areaerrorbar_wosim(t,norm_sigH,norm_sigH_std,mycolors(3,:))% error bar for real data
hold on
plot_all_sim_2(t,norm_H_conc_t(1:end-1)',norm_H_conc_t_std(1:end-1)',mycolors(4,:)) % color patch for predicted data
% legend('Experimental \sigma^{A}','','Predicted \sigma^{A}','Experimental \sigma^{H}','','Predicted \sigma^{H}','location','northwest')
legend('Phy-spank (\sigma^{A})','','spo0A Pv','PcitG (\sigma^{H})','','spo0A Ps','location','northwest','FontSize',22)

xlabel('Time (h)')
ylabel('Normalized Promoter Activity')
set(gca,'FontSize',28)

% xlim([2 9])
% ylim([0 5600])
% legend('WT')
f = gcf;
exportgraphics(f,'SigmaDynamics_allpromtoers_norm.tif','Resolution',600)

