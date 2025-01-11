clear
clc
best_fit_model_prediction

%%
close all
clear
clc
load('sim_data_original.mat')
load('new_promoter_activity_single.mat');
% load('random_dataset.mat');

real_data_Ps = new_real_data_Ps(1:end-1,:); % exclude the real WT, instead WT is Ps only, no mutations
% real_data_Ps = Ps_rand;
real_data_Ps_std = new_Ps_only_std(1:end-1,:);

real_data_Pv = new_real_data_Pv(1:end-1,:);
% real_data_Pv = Pv_rand;
real_data_Pv_std = new_Pv_only_std(1:end-1,:);

load('new_PsPv_promoter_activity.mat')
real_data_PsPv = PvPs_avg_new(1:end-1,:);
real_data_Psv_std = PvPs_std_new(1:end-1,:);

TF_conc_t = get_aps(2:9);
%% see if two-promoter data is fit well
load('best_fit_model_prediction_colors.mat')
rgbColor = hex2rgb(hexColor);

% sim_all = sim_data_v+sim_data_s;
sim_all = sim_data_sv;

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim(fliplr(real_data_PsPv(:,[8,1:4])),fliplr(real_data_Psv_std(:,[8,1:4])),flipud(rgbColor([8,1:4],:)))
hold on
plot_all_sim(fliplr(sim_all(:,[8,1:4])),flipud(rgbColor([8,1:4],:)))
% legend('','Real WT','Sim WT','','Real 1*','Sim 1*','','Real 2*','Sim 2*','','Real 3*','Sim 3*','Location','northwest','FontSize',16,'NumColumns',2)
% legend('','None','','','1*','','','2*','','','3*','','Location','northwest','FontSize',22,'NumColumns',2)
% legend('','','','','','','','','None','1*','2*','3*','Location','northwest','FontSize',22,'NumColumns',2)
% with patch only
% legend('None','1*','2*','3*','Location','northwest','FontSize',22,'NumColumns',2)
% with errorbar only
[~, hobj, ~, ~] = legend('','','','','','12','13','23','123','none','Location','northwest','FontSize',22,'NumColumns',3);
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.87);


% title('PvPs together, Single Mutation')
xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
xlim([2 9])
ylim([0 1400])
set(gca,'FontSize',22)
f = gcf;
exportgraphics(f,'Pure_kinetic_PsPv_sing.tif','Resolution',600)

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim(fliplr(real_data_PsPv(:,[8,1,5:7])),fliplr(real_data_Psv_std(:,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
hold on
plot_all_sim(fliplr(sim_all(:,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
% title('PvPs together, Multiple Mutations')
% legend('','None','','','12*','','','13*','','','23*','','','123*','','Location','northwest','FontSize',22,'NumColumns',2)
% legend('','','','','','','','','','','None','12*','13*','23*','123*','Location','northwest','FontSize',22,'NumColumns',2)
% with patch only
% legend('None','12*','13*','23*','Location','northwest','FontSize',22,'NumColumns',2)
% with errorbar only
[~,hobj,~,~] = legend('','','','','','1','2','3','123','none','Location','northwest','FontSize',22,'NumColumns',3);
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.87);
xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
xlim([2 9])
ylim([0 1400])
set(gca,'FontSize',22)

% saveas(gcf,'real_data_ratio_vs_time.fig')
f = gcf;
exportgraphics(f,'Pure_kinetic_PsPv_muti.tif','Resolution',600)

fc = (sim_all(end,:)-sim_all(end,end))./sim_all(end,end);

%% Ps real and simulation on the same graph, unnormalized
% unnorm,type 1234
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim(fliplr(real_data_Ps(:,[8,1:4])),fliplr(real_data_Ps_std(:,[8,1:4])),flipud(rgbColor([8,1:4],:)))
hold on
plot_all_sim(fliplr(sim_data_s(:,[8,1:4])),flipud(rgbColor([8,1:4],:)))
% title('Ps only, Single Mutations')
% legend('','Real WT','Sim WT','','Real 1*','Sim 1*','','Real 2*','Sim 2*','','Real 3*','Sim 3*','Location','northwest','NumColumns',2,'FontSize',16)
% legend('','None','','','1*','','','2*','','','3*','','Location','northwest','FontSize',22,'NumColumns',2)
% legend('','','','','','','','','None','1*','2*','3*','Location','northwest','FontSize',22,'NumColumns',2)
% legend('','','','','','','','','None','1*','2*','3*','Location','northwest','FontSize',22,'NumColumns',2)
% legend('None','1*','2*','3*','Location','northwest','FontSize',22,'NumColumns',2)
[~, hobj, ~, ~] = legend('','','','','','12','13','23','123','none','Location','northwest','FontSize',22,'NumColumns',3,'Orientation','vertical');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.87);
xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
ylim([0 1400])
xlim([2 9])
set(gca,'FontSize',22)
f = gcf;
exportgraphics(f,'Pure_kinetic_Ps_sing.tif','Resolution',600)

% unnorm,type 15678
% subplot(2,2,2)
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim(fliplr(real_data_Ps(:,[8,1,5:7])),fliplr(real_data_Ps_std(:,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
hold on
plot_all_sim(fliplr(sim_data_s(:,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
% title('Ps only, Multiple Mutations')
% legend('','Real WT','Sim WT','','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*','Location','northwest','NumColumns',2,'FontSize',12)
% legend('','None','','','12*','','','13*','','','23*','','','123*','','Location','northwest','FontSize',22,'NumColumns',3)
% legend('','','','','','','','','','','None','12*','13*','23*','123*','Location','northwest','FontSize',22,'NumColumns',2)
% legend('None','12*','13*','23*','Location','northwest','FontSize',22,'NumColumns',2)
[~,hobj,~,~] = legend('','','','','','1','2','3','123','none','Location','northwest','FontSize',22,'NumColumns',3,'Orientation','vertical');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.87);
xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
ylim([0 1400])
xlim([2 9])
set(gca,'FontSize',22)

f = gcf;
exportgraphics(f,'Pure_kinetic_Ps_muti.tif','Resolution',600)
fc = (sim_data_s(end,:)-sim_data_s(end,end))./sim_data_s(end,end);

%%
% % plot just wt and 123* to see why there is a big difference etween
% % datacollapsed model and datacollapsed real data
% figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
% plot_areaerrorbar_wosim(fliplr(real_data_Ps(:,[8,1])),fliplr(real_data_Ps_std(:,[8,1])),flipud(rgbColor([8,1],:)))
% hold on
% plot_all_sim(fliplr(sim_data_s(:,[8,1])),flipud(rgbColor([8,1],:)))
% % title('Ps only, Multiple Mutations')
% % legend('','Real WT','Sim WT','','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*','Location','northwest','NumColumns',2,'FontSize',12)
% % legend('','None','','','12*','','','13*','','','23*','','','123*','','Location','northwest','FontSize',22,'NumColumns',3)
% % legend('','','','','','','','','','','None','12*','13*','23*','123*','Location','northwest','FontSize',22,'NumColumns',2)
% % legend('None','12*','13*','23*','Location','northwest','FontSize',22,'NumColumns',2)
% [~,hobj,~,~] = legend('','','123','none','Location','northwest','FontSize',22,'NumColumns',3,'Orientation','vertical');
% hl = findobj(hobj,'type','line');
% set(hl,'LineWidth',2.87);
% xlabel('Time (h)')
% ylabel('v^{effective} (a.u.)')
% ylim([0 650])
% xlim([2 9])
% set(gca,'FontSize',22)

%% Pv real and simulation on the same graph, unnormalized
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
% subplot(2,1,1)
plot_areaerrorbar_wosim(fliplr(real_data_Pv(:,[8,1:4])),fliplr(real_data_Pv_std(:,[8,1:4])),flipud(rgbColor([8,1:4],:)));
hold on
plot_all_sim(fliplr(sim_data_v(:,[8,1:4])),flipud(rgbColor([8,1:4],:)))
% title('Pv only, Single Mutations')
% legend('','Real WT','Sim WT','','Real 1*','Sim 1*','','Real 2*','Sim 2*','','Real 3*','Sim 3*','Location','northwest','NumColumns',2,'FontSize',16)
% legend('','None','','','1*','','','2*','','','3*','','Location','northwest','FontSize',22,'NumColumns',3)
% legend('','','','','','','','','None','1*','2*','3*','Location','northwest','FontSize',22,'NumColumns',2)
% legend('None','1*','2*','3*','Location','northwest','FontSize',22,'NumColumns',2)
[~, hobj, ~, ~] = legend('','','','','','12','13','23','123','none','Location','northwest','FontSize',22,'NumColumns',3,'Orientation','vertical');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.87);

xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
ylim([0 1400])
xlim([2 9])
set(gca,'FontSize',22)
f = gcf;
exportgraphics(f,'Pure_kinetic_Pv_sing.tif','Resolution',600)

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim(fliplr(real_data_Pv(:,[8,1,5:7])),fliplr(real_data_Pv_std(:,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
hold on
plot_all_sim(fliplr(sim_data_v(:,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
% title('Pv only, Multiple Mutations')
% legend('','Real WT','Sim WT','','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*','Location','northwest','NumColumns',2,'FontSize',12)
% legend('','None','','','12*','','','13*','','','23*','','','123*','','Location','northwest','FontSize',22,'NumColumns',3)
% legend('','','','','','','','','','','None','12*','13*','23*','123*','Location','northwest','FontSize',22,'NumColumns',2)
% legend('None','12*','13*','23*','Location','northwest','FontSize',22,'NumColumns',2)
[~,hobj,~,~] = legend('','','','','','1','2','3','123','none','Location','northwest','FontSize',22,'NumColumns',3,'Orientation','vertical');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.87);
xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
xlim([2 9])
ylim([0 1400])
set(gca,'FontSize',22)
f = gcf;
exportgraphics(f,'Pure_kinetic_Pv_muti.tif','Resolution',600)
fc = (sim_data_v(end,:)-sim_data_v(end,end))./sim_data_v(end,end);