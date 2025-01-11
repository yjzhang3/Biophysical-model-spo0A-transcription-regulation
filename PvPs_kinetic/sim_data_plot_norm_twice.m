% normalize the data twice, by first time point, and by none data (data
% collapse)
%%%% load real data
clear
best_fit_model_prediction
load('best_fit_model_prediction_colors.mat')
load('sim_data_original.mat')
close all
%% Normalize real data
nn_real_s = sim_data_s./sim_data_s(:,end);
nn_real_s_final = nn_real_s./nn_real_s(8,:);

nn_real_v = sim_data_v./sim_data_v(:,end);
nn_real_v_final = nn_real_v./nn_real_v(8,:);

nn_real_sv = sim_data_sv./sim_data_sv(:,end);
nn_real_sv_final = nn_real_sv./nn_real_sv(8,:);

rgbColor = hex2rgb(hexColor);

%% plotting norm real data Ps
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
rgbColor = hex2rgb(hexColor);

% subplot(2,1,1)
plot_all_sim(fliplr(nn_real_s_final(1:end,[8,1:4])),flipud(rgbColor([8,1:4],:)))
% legend('','None','','1*','','2*','','3*','NumColumns',2,'Location','Northwest')
legend('12','13','23','123','none','Location','northeast','FontSize',20,'NumColumns',3,'Orientation','horizontal')

% title('Ps Strain, Single Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2])
xlim([2 9])
set(gca,'FontSize',22)

% saveas(gcf,'Psonly_real_sing_norm.fig')
f = gcf;
exportgraphics(f,'Psonly_sim_sing_norm_twice.tif','Resolution',600)

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_all_sim(fliplr(nn_real_s_final(1:end,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
% legend('','None','','12*','','13*','','23*','','123*','NumColumns',2,'Location','Northwest')
legend('1','2','3','123','none','Location','northeast','FontSize',20,'NumColumns',3,'Orientation','horizontal')

% title('Ps Strain, Multiple Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2])
xlim([2 9])
set(gca,'FontSize',22)


% saveas(gcf,'Psonly_real_muti_norm.fig')
% saveas(gcf,'Psonly_real_muti.tif')
f = gcf;
exportgraphics(f,'Psonly_sim_muti_norm_twice.tif','Resolution',600)
%% plotting norm real data Pv
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
% hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#39F9D9';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

% subplot(2,1,1)
plot_all_sim(fliplr(nn_real_v_final(1:end,[8,1:4])),flipud(rgbColor([8,1:4],:)))
legend('12','13','23','123','none','Location','northeast','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('Pv Strain, Single Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2])
xlim([2 9])
set(gca,'FontSize',22)

% saveas(gcf,'Pvonly_real_sing_norm.fig')
f = gcf;
exportgraphics(f,'Pvonly_sim_sing_norm_twice.tif','Resolution',600)

% subplot(2,1,2)
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_all_sim(fliplr(nn_real_v_final(1:end,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
legend('1','2','3','123','none','Location','northeast','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('Pv Strain, Multiple Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2])
xlim([2 9])
set(gca,'FontSize',22)


% saveas(gcf,'Pvonly_real_muti_norm.fig')
f = gcf;
exportgraphics(f,'Pvonly_sim_muti_norm_twice.tif','Resolution',600)

%% plotting norm real data PsPv
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
% hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#39F9D9';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

% subplot(2,1,1)
plot_all_sim(fliplr(nn_real_sv_final(1:end,[8,1:4])),flipud(rgbColor([8,1:4],:)))
legend('12','13','23','123','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('PvPs Strain, Single Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2])
xlim([2 9])
set(gca,'FontSize',22)

% saveas(gcf,'Psv_real_sing_norm_twice.fig')
f = gcf;
% exportgraphics(f,'Psv_real_sing.tif','Resolution',600)
exportgraphics(f,'Psv_sim_sing_norm_twice.tif','Resolution',600)

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
plot_all_sim(fliplr(nn_real_sv_final(1:end,[8,1,5:7])),flipud(rgbColor([8,1,5:7],:)))
% legend('','None','','12*','','13*','','23*','','123*','NumColumns',2,'Location','Northwest')
legend('1','2','3','123','none','Location','northwest','FontSize',20,'NumColumns',3,'Orientation','horizontal')
% title('PvPs Strain, Multiple Mutations')
xlabel('Time (h)')
ylabel('Normalized promoter activity')
ylim([0 2])
xlim([2 9])
set(gca,'FontSize',22)


% saveas(gcf,'Psv_real_muti_norm_twice.fig')
f = gcf;
exportgraphics(f,'Psv_sim_muti_norm_twice.tif','Resolution',600)

% %% compute maximum fold change in unnormalized data and normalized data 
% s_max_fold_change_per_strain_unnorm = max(real_data_Ps,[],1)./min(real_data_Ps,[],1);
% s_max_fold_change_per_strain_norm = max(nn_real_s_final,[],1)./min(nn_real_s_final,[],1);
% 
% v_max_fold_change_per_strain_unnorm = max(real_data_Pv,[],1)./min(real_data_Pv,[],1);
% v_max_fold_change_per_strain_norm = max(nn_real_v_final,[],1)./min(nn_real_v_final,[],1);
% 
% sv_max_fold_change_per_strain_unnorm = max(real_data_PsPv,[],1)./min(real_data_PsPv,[],1);
% sv_max_fold_change_per_strain_norm = max(nn_real_sv_final,[],1)./min(nn_real_sv_final,[],1);
