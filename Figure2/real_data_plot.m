%% load real data
clear
load('best_fit_model_prediction_colors.mat')
rgbColor = hex2rgb(hexColor);

load('new_promoter_activity_single.mat');
% load('random_dataset.mat');

real_data_Ps = new_real_data_Ps; % exclude the real WT, instead WT is Ps only, no mutations
% real_data_Ps = Ps_rand;
real_data_Ps_std = new_Ps_only_std;

real_data_Pv = new_real_data_Pv;
% real_data_Pv = Pv_rand;
real_data_Pv_std = new_Pv_only_std;

load('new_PsPv_promoter_activity.mat')
real_data_PsPv = PvPs_avg_new;
real_data_PsPv_std = PvPs_std_new;

%% fold change between two strains at t = 2h
r = real_data_Pv(1,1)/real_data_Ps(1,1);

%% fold change for each strain
v_change = real_data_Pv(1,end)/real_data_Pv(1,1);
s_change = real_data_Ps(1,end)/real_data_Ps(1,1);

%% real data exaggerate
% figure('position',[188,100,1177/2,510]); % for additivity
% % errorbar(real_sum(1:end-1,:),real_sum_std(1:end-1,:),'Marker','o','MarkerSize',4,'LineWidth',2)
% % ax = gca;
% % ax.ColorOrder = rgbColor;
% % ax.FontSize = 20;
% hexColor=['#000000';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
% rgbColor = hex2rgb(hexColor);
% 
% 
% plot_areaerrorbar_wosim_2(real_sum(1:end-1,:),real_sum_std(1:end-1,:),rgbColor)
% legend('','None','','1*','','2*','','3*','','12*','','13*','','23*','','123*','NumColumns',4)
% % title('Additivity in Real Data')
% xlabel('Time (h)')
% ylabel('Additivity Level')
% ylim([0 2.3])
% xlim([2 9])
% set(gca,'FontSize',22)
% 
% saveas(gcf,'real_data_ratio_vs_time.fig')
% f = gcf;
% exportgraphics(f,'real_data_ratio_vs_time.tif','Resolution',600)


%% plotting real data Ps
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
% hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#39F9D9';'#ffa300';'#dc0ab4';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

% subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(real_data_Ps(1:end-1,[8,1:4])),fliplr(real_data_Ps_std(1:end-1,[8,1:4])),flipud(rgbColor([8,1:4],:)))
% legend('','None','','1*','','2*','','3*','NumColumns',2,'Location','Northwest')
legend('','12','','13','','23','','123','','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')

% title('Ps Strain, Single Mutations')
xlabel('Time (h)')
ylabel('Promoter activity (Miller units)')
ylim([0 1400])
xlim([2 9])
set(gca,'FontSize',22)

saveas(gcf,'Psonly_real_sing.fig')
f = gcf;
exportgraphics(f,'Psonly_real_sing.tif','Resolution',300)

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(real_data_Ps(1:end-1,[8,1,5:7])),fliplr(real_data_Ps_std(1:end-1,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
% legend('','None','','12*','','13*','','23*','','123*','NumColumns',2,'Location','Northwest')
legend('','1','','2','','3','','123','','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')

% title('Ps Strain, Multiple Mutations')
xlabel('Time (h)')
ylabel('Promoter activity (Miller units)')
ylim([0 1400])
xlim([2 9])
set(gca,'FontSize',22)


saveas(gcf,'Psonly_real_muti.fig')
% saveas(gcf,'Psonly_real_muti.tif')
f = gcf;
exportgraphics(f,'Psonly_real_muti.tif','Resolution',300)
%% plotting real data Pv
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
% hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#39F9D9';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

% subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(real_data_Pv(1:end-1,[8,1:4])),fliplr(real_data_Pv_std(1:end-1,[8,1:4])),flipud(rgbColor([8,1:4],:)))
legend('','12','','13','','23','','123','','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('Pv Strain, Single Mutations')
xlabel('Time (h)')
ylabel('Promoter activity (Miller units)')
ylim([0 1400])
xlim([2 9])
set(gca,'FontSize',22)

saveas(gcf,'Pvonly_real_sing.fig')
f = gcf;
exportgraphics(f,'Pvonly_real_sing.tif','Resolution',600)

% subplot(2,1,2)
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(real_data_Pv(1:end-1,[8,1,5:7])),fliplr(real_data_Pv_std(1:end-1,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
legend('','1','','2','','3','','123','','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('Pv Strain, Multiple Mutations')
xlabel('Time (h)')
ylabel('Promoter activity (Miller units)')
ylim([0 1400])
xlim([2 9])
set(gca,'FontSize',22)


saveas(gcf,'Pvonly_real_muti.fig')
f = gcf;
exportgraphics(f,'Pvonly_real_muti.tif','Resolution',600)

%% plotting real data PsPv
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
% hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#39F9D9';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

% subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(real_data_PsPv(1:end-1,[8,1:4])),fliplr(real_data_PsPv_std(1:end-1,[8,1:4])),flipud(rgbColor([8,1:4],:)))
legend('','12','','13','','23','','123','','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('PvPs Strain, Single Mutations')
xlabel('Time (h)')
ylabel('Promoter activity (Miller units)')
ylim([0 1400])
xlim([2 9])
set(gca,'FontSize',22)

saveas(gcf,'Psv_real_sing.fig')
f = gcf;
% exportgraphics(f,'Psv_real_sing.tif','Resolution',600)
exportgraphics(f,'Psv_real_sing.tif','Resolution',600)

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(real_data_PsPv(1:end-1,[8,1,5:7])),fliplr(real_data_PsPv_std(1:end-1,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
% legend('','None','','12*','','13*','','23*','','123*','NumColumns',2,'Location','Northwest')
legend('','1','','2','','3','','123','','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('PvPs Strain, Multiple Mutations')
xlabel('Time (h)')
ylabel('Promoter activity (Miller units)')
ylim([0 1400])
xlim([2 9])
set(gca,'FontSize',22)


saveas(gcf,'Psv_real_muti.fig')
f = gcf;
exportgraphics(f,'Psv_real_muti.tif','Resolution',600)
