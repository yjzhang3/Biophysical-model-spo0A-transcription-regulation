% normalize the data twice, by last time point, and by none data (data
% collapse)
%% load real data
clear
clc
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

real_sum = (real_data_Pv+real_data_Ps)./PvPs_avg_new;
real_sum_std = (real_data_Ps_std + real_data_Pv_std)./PvPs_avg_new;

%% Normalize real data
tp = 8; % by which time point do we normalize
nn_real_s = real_data_Ps./real_data_Ps(:,end);
nn_real_s_final = nn_real_s./nn_real_s(tp,:);
nn_real_s_std = real_data_Ps_std./real_data_Ps(:,end);
nn_real_s_std_final = nn_real_s_std./nn_real_s(tp,:);

nn_real_v = real_data_Pv./real_data_Pv(:,end);
nn_real_v_final = nn_real_v./nn_real_v(tp,:);
nn_real_v_std = real_data_Pv_std./real_data_Pv(:,end);
nn_real_v_std_final = nn_real_v_std./nn_real_v(tp,:);

nn_real_sv = real_data_PsPv./real_data_PsPv(:,end);
nn_real_sv_final = nn_real_sv./nn_real_sv(tp,:);
nn_real_sv_std = real_data_PsPv_std./real_data_PsPv(:,end);
nn_real_sv_std_final = nn_real_sv_std./nn_real_sv(tp,:);
hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#39F9D9';'#ffa300';'#dc0ab4';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

%% plotting norm real data Ps
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
rgbColor = hex2rgb(hexColor);

% subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(nn_real_s_final(1:end-1,[1:4])),fliplr(nn_real_s_std_final(1:end-1,[1:4])),flipud(rgbColor([1:4],:)))
% legend('','None','','1*','','2*','','3*','NumColumns',2,'Location','Northwest')
legend('','12','','13','','23','','123','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')

% title('Ps Strain, Single Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2.2])
xlim([2 9])
set(gca,'FontSize',22)

% saveas(gcf,'Psonly_real_sing_norm.fig')
f = gcf;
exportgraphics(f,'Psonly_real_sing_norm_twice.tif','Resolution',600)

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(nn_real_s_final(1:end-1,[8,5:7])),fliplr(nn_real_s_std_final(1:end-1,[8,5:7])),flipud(rgbColor([8,5:7],:)))
% legend('','None','','12*','','13*','','23*','','123*','NumColumns',2,'Location','Northwest')
legend('','1','','2','','3','','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')

% title('Ps Strain, Multiple Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2.2])
xlim([2 9])
set(gca,'FontSize',22)


% saveas(gcf,'Psonly_real_muti_norm.fig')
% saveas(gcf,'Psonly_real_muti.tif')
f = gcf;
exportgraphics(f,'Psonly_real_muti_norm_twice.tif','Resolution',600)
%% plotting norm real data Pv
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
% hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#39F9D9';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

% subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(nn_real_v_final(1:end-1,[1:4])),fliplr(nn_real_v_std_final(1:end-1,[1:4])),flipud(rgbColor([1:4],:)))
legend('','12','','13','','23','','123','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('Pv Strain, Single Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2.2])
xlim([2 9])
set(gca,'FontSize',22)

% saveas(gcf,'Pvonly_real_sing_norm.fig')
f = gcf;
exportgraphics(f,'Pvonly_real_sing_norm_twice.tif','Resolution',600)

% subplot(2,1,2)
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(nn_real_v_final(1:end-1,[8,5:7])),fliplr(nn_real_v_std_final(1:end-1,[8,5:7])),flipud(rgbColor([8,5:7],:)))
legend('','1','','2','','3','','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('Pv Strain, Multiple Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2.2])
xlim([2 9])
set(gca,'FontSize',22)


% saveas(gcf,'Pvonly_real_muti_norm.fig')
f = gcf;
exportgraphics(f,'Pvonly_real_muti_norm_twice.tif','Resolution',600)

%% plotting norm real data PsPv
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
% hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#39F9D9';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

% subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(nn_real_sv_final(1:end-1,[1:4])),fliplr(nn_real_sv_std_final(1:end-1,[1:4])),flipud(rgbColor([1:4],:)))
legend('','12','','13','','23','','123','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('PvPs Strain, Single Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2.2])
xlim([2 9])
set(gca,'FontSize',22)

% saveas(gcf,'Psv_real_sing_norm_twice.fig')
f = gcf;
% exportgraphics(f,'Psv_real_sing.tif','Resolution',600)
exportgraphics(f,'Psv_real_sing_norm_twice.tif','Resolution',600)

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_areaerrorbar_wosim_2(fliplr(nn_real_sv_final(1:end-1,[8,5:7])),fliplr(nn_real_sv_std_final(1:end-1,[8,5:7])),flipud(rgbColor([8,5:7],:)))
% legend('','None','','12*','','13*','','23*','','123*','NumColumns',2,'Location','Northwest')
legend('','1','','2','','3','','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('PvPs Strain, Multiple Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2.2])
xlim([2 9])
set(gca,'FontSize',22)


% saveas(gcf,'Psv_real_muti_norm_twice.fig')
f = gcf;
exportgraphics(f,'Psv_real_muti_norm_twice.tif','Resolution',600)

%% compute maximum fold change in unnormalized data and normalized data 
s_max_fold_change_per_strain_unnorm = max(real_data_Ps,[],1)./min(real_data_Ps,[],1);
s_max_fold_change_per_strain_norm = max(nn_real_s_final,[],1)./min(nn_real_s_final,[],1);

v_max_fold_change_per_strain_unnorm = max(real_data_Pv,[],1)./min(real_data_Pv,[],1);
v_max_fold_change_per_strain_norm = max(nn_real_v_final,[],1)./min(nn_real_v_final,[],1);

sv_max_fold_change_per_strain_unnorm = max(real_data_PsPv,[],1)./min(real_data_PsPv,[],1);
sv_max_fold_change_per_strain_norm = max(nn_real_sv_final,[],1)./min(nn_real_sv_final,[],1);

%% plot real data additivity level
figure('position',[188,100,1300,550]); % for additivity
% errorbar(real_sum(1:end-1,:),real_sum_std(1:end-1,:),'Marker','o','MarkerSize',4,'LineWidth',2)
% ax = gca;
% ax.ColorOrder = rgbColor;
% ax.FontSize = 20;
hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#39F9D9';'#ffa300';'#dc0ab4';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);


plot_areaerrorbar_wosim_2(fliplr(real_sum(1:end-1,[8,1,2:7])),fliplr(real_sum_std(1:end-1,[8,1,2:7])),flipud(rgbColor([8,1,2:7],:)))
legend('','1','','2','','3','','12','','13','','23','','123','','None','NumColumns',4,'FontSize',20,'location','Northwest','Orientation','horizontal')
% title('Additivity in Real Data')
xlabel('Time (h)')
ylabel('Ratio')
ylim([0 2.3])
xlim([2 9])
set(gca,'FontSize',22)

saveas(gcf,'real_data_ratio_vs_time.fig')
f = gcf;
exportgraphics(f,'real_data_ratio_vs_time.tif','Resolution',600)
