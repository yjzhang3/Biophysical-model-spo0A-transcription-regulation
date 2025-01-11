clear 
clc
close all
fn = 'purethermo_Sept_2024_extrarun_highera.txt';
A = readmatrix(fn);
[M,I] = min(A(:,end));
pars = A(I,:);
% pars = A(1,:);
pars = pars(1:end-1);
% pars = [0.037174,1.000000,1.000000,0.101697,0.006738,1.052096,0.641727,20.085177,0.613681,0.739486,1.412397,0.778704,1.000000,8.148579,1.068196,0.654056,1.904700,0.814109,1.000000,0.110578,1.757030,1.000000,1.000000,1,1,1,1,1];

%%
strain_type = 1:8;
% group_array_s.g1 = 1:16;
group_array_s.g1 = 1:8;
group_array_v.g1 = 1:8;
group_array_sv.g1 = 1:8;

load('new_promoter_activity_single.mat');
load('new_PsPv_promoter_activity.mat');

real_data_Ps = new_real_data_Ps(1:end-1,strain_type); % exclude the real WT, instead WT is Ps only, no mutations
% real_data_Ps = Ps_rand;
real_data_Ps_std = new_Ps_only_std(1:end-1,strain_type);

real_data_Pv = new_real_data_Pv(1:end-1,strain_type);
% real_data_Pv = Pv_rand;
real_data_Pv_std = new_Pv_only_std(1:end-1,strain_type);

real_data_PsPv = PvPs_avg_new(1:end-1,strain_type);
real_data_PsPv_std = PvPs_std_new(1:end-1,strain_type);

TF_conc_t = get_aps(2:9);

% nbd = 6;
nbd = 5;
mut_mat_s = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_v = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_sv = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_s = mut_mat_s(strain_type,:);
mut_mat_v = mut_mat_v(strain_type,:);
mut_mat_sv = mut_mat_sv(strain_type,:);

%%
 [final_diff,sim_data_Ps,sim_data_Pv,sim_data_PsPv,vmax_s,vmax_v,vmax_sv,...
     vmax_array_s,vmax_array_v,vmax_array_sv,...
    energyi_s,energyi_v,energyi_sv,H_conc_t,A_conc_t] = objective_function(...
    nbd,pars,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv);


norm_energyi_G = -log(energyi_sv);

load('best_fit_model_prediction_colors.mat')
rgbColor_raw = hex2rgb(hexColor);
rgbColor = rgbColor_raw(strain_type,:);

%% plot probability of each configuration contributing to a strain over time

strain_name_s = {'Ps WT','Ps 1*','Ps 2*','Ps 3*','Ps 12*','Ps 13*','Ps 23*','Ps 123*'};
strain_name_v = {'Pv WT','Pv 1*','Pv 2*','Pv 3*','Pv 12*','Pv 13*','Pv 23*','Pv 123*'};
strain_name_sv = {'PsPv WT','PsPv 1*','PsPv 2*','PsPv 3*','PsPv 12*','PsPv 13*','PsPv 23*','PsPv 123*'};


%% test additivity
figure("position",[680,678,560,420]); % for additivity
ratio = (sim_data_Pv+sim_data_Ps)./sim_data_PsPv;
plot(2:9,ratio,'LineWidth',2)
ylim([0 2])
xlim([2 9])
hold on
legend('WT','1*','2*','3*','12*','13*','23*','123*','NumColumns',2)
xlabel('Time (hours)')
ylabel('Ratio')

load('best_fit_model_prediction_colors.mat')
rgbColor = hex2rgb(hexColor);
set(gca,'FontSize',16,'ColorOrder', rgbColor)

%% Ps real and simulation on the same graph, unnormalized and normalized
% hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
% rgbColor = hex2rgb(hexColor);

ind_single = find(strain_type<5);
ind_mul = find(strain_type>4);

% unnorm,type 1234
figure('position',[188,174,1177,510]); % for strain fits, single row
% subplot(2,2,1)+`
subplot(1,2,1)
plot_areaerrorbar_wsim(real_data_Ps(:,ind_single),real_data_Ps_std(:,ind_single),sim_data_Ps(:,ind_single),rgbColor(ind_single,:));
title('Ps-only, Single Mutations')
A = strain_type(ind_single)*3-2;
B = A+1;
C = B+1;
AB = reshape([A;B;C], size(A,1), []);
led = {'','Real WT','Sim WT','','Real 1*','Sim 1*','','Real 2*','Sim 2*','','Real 3*','Sim 3*'};
legend(led(AB),'Location','northwest','NumColumns',2,'Location','Northwest','FontSize',16)
xlabel('Time (hours)')
ylabel('Transcription Rate')
ylim([0 1200])
xlim([2 9])
set(gca,'FontSize',20)

% unnorm,type 15678 % wt is also plotted!
% subplot(2,2,2)
% subplot(1,2,2)
% plot_areaerrorbar_wsim(real_data_Ps(:,[1,5:8]),real_data_Ps_std(:,[1,5:8]),sim_data_Ps(:,[1,5:8]),rgbColor([1,5:8],:))
% title('Ps-only, Multiple Mutations')
% legend('','Real WT','Sim WT','','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*','Location','northwest','NumColumns',2,'Location','Northwest','FontSize',12)
% xlabel('Time (hours)')
% ylabel('Transcription Rate')
% ylim([0 1200])
% xlim([2 9])
% set(gca,'FontSize',20)

subplot(1,2,2)
plot_areaerrorbar_wsim(real_data_Ps(:,ind_mul),real_data_Ps_std(:,ind_mul),sim_data_Ps(:,ind_mul),rgbColor(ind_mul,:))
title('Ps-only, Multiple Mutations')
A = strain_type(ind_mul)*3-2;
B = A+1;
C = B+1;
AB = reshape([A;B;C], size(A,1), []);
led={'','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*'};
legend(led(AB-12),'Location','northwest','NumColumns',2,'Location','Northwest','FontSize',12);
xlabel('Time (hours)')
ylabel('Transcription Rate')
ylim([0 1200])
xlim([2 9])
set(gca,'FontSize',20)
saveas(gcf,'Pure_thermo_quaE_Ps.jpg')
saveas(gcf,'Pure_thermo_quaE_Ps.fig')

%% Pv real and simulation on the same graph unnormalized and normalized
figure('position',[188,174,1177,510]); % for strain fits, single row

ind_single = find(strain_type<5);
ind_mul = find(strain_type>4);

% subplot(2,2,1)
subplot(1,2,1)
plot_areaerrorbar_wsim(real_data_Pv(:,ind_single),real_data_Pv_std(:,ind_single),sim_data_Pv(:,ind_single),rgbColor(ind_single,:));
A = strain_type(ind_single)*3-2;
B = A+1;
C = B+1;
AB = reshape([A;B;C], size(A,1), []);
title('Pv-only, Single Mutations')
led = {'','Real WT','Sim WT','','Real 1*','Sim 1*','','Real 2*','Sim 2*','','Real 3*','Sim 3*','Location','northwest'};
legend(led(AB),'NumColumns',2,'Location','Northwest','FontSize',16)
xlabel('Time (hours)')
ylabel('Transcription Rate')
ylim([0 500])
xlim([2 9])
set(gca,'FontSize',20)

% subplot(2,2,2)
% plot_areaerrorbar_wsim(real_data_Pv(:,[1,5:8]),real_data_Pv_std(:,[1,5:8]),sim_data_Pv(:,[1,5:8]),rgbColor([1,5:8],:))
% title('Pv-only, Multiple Mutations')
% legend('','Real WT','Sim WT','','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*','Location','northwest','NumColumns',2,'Location','Northwest','FontSize',12)
% xlabel('Time (hours)')
% ylabel('Transcription Rate')
% ylim([0 500])
% xlim([2 9])
% set(gca,'FontSize',20)

subplot(1,2,2)
plot_areaerrorbar_wsim(real_data_Pv(:,ind_mul),real_data_Pv_std(:,ind_mul),sim_data_Pv(:,ind_mul),rgbColor(ind_mul,:))
A = strain_type(ind_mul)*3-2;
B = A+1;
C = B+1;
AB = reshape([A;B;C], size(A,1), []);
title('Pv-only, Multiple Mutations')
led = {'','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*'};
legend(led(AB-12),'Location','northwest','NumColumns',2,'Location','Northwest','FontSize',12)
xlabel('Time (hours)')
ylabel('Transcription Rate')
ylim([0 500])
xlim([2 9])
set(gca,'FontSize',20)
saveas(gcf,'Pure_thermo_quaE_Pv.jpg')
saveas(gcf,'Pure_thermo_quaE_Pv.fig')

%% see two promoter model fits to pspv data well
% hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
% rgbColor = hex2rgb(hexColor);
% load('new_PsPv_promoter_activity.mat')

sim_all = sim_data_PsPv;

ind_single = find(strain_type<5);
ind_mul = find(strain_type>4);

figure('position',[188,174,1177,510]); % for strain fits, single row
subplot(1,2,1)
plot_areaerrorbar_wsim(real_data_PsPv(:,ind_single),real_data_PsPv_std(:,ind_single),sim_all(:,ind_single),rgbColor(ind_single,:))
A = strain_type(ind_single)*3-2;
B = A+1;
C = B+1;
AB = reshape([A;B;C], size(A,1), []);
led = {'','Real WT','Sim WT','','Real 1*','Sim 1*','','Real 2*','Sim 2*','','Real 3*','Sim 3*'};
legend(led(AB),'NumColumns',2,'Location','Northwest','FontSize',16)
title('PsPv Together, Single Mutation','FontSize',20)
xlabel('Time (hours)')
ylabel('Transcription Rate')
xlim([2 9])
ylim([0 1400])
set(gca,'FontSize',20)

subplot(1,2,2)
% plot_areaerrorbar_wsim(PvPs_avg_new(1:end-1,[1,5:8]),PvPs_std_new(1:end-1,[1,5:8]),sim_all(:,[1,5:8]),rgbColor([1,5:8],:))
% title('PsPv Together, Multiple Mutations','FontSize',20)
% legend('','Real WT','Sim WT','','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*','NumColumns',2,'Location','Northwest','FontSize',12)
% xlabel('Time (hours)')
% ylabel('Transcription Rate')
% xlim([2 9])
% ylim([0 1400])
% set(gca,'FontSize',20)

plot_areaerrorbar_wsim(real_data_PsPv(:,ind_mul),real_data_PsPv_std(:,ind_mul),sim_all(:,ind_mul),rgbColor(ind_mul,:))
title('PsPv Together, Multiple Mutations','FontSize',20)
A = strain_type(ind_mul)*3-2;
B = A+1;
C = B+1;
AB = reshape([A;B;C], size(A,1), []);
led = {'','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*'};
legend(led(AB-12),'NumColumns',2,'Location','Northwest','FontSize',12)
xlabel('Time (hours)')
ylabel('Transcription Rate')
xlim([2 9])
ylim([0 1400])
set(gca,'FontSize',20)
saveas(gcf,'Pure_thermo_quaE_PsPv.jpg')
saveas(gcf,'Pure_thermo_quaE_PsPv.fig')

%% obtain vmax and energy for multiple runs
A = readmatrix(fn);
energyi_sv_multiple_run = zeros(length(A(:,1)),length(energyi_sv));
energyi_s_multiple_run = zeros(length(A(:,1)),length(energyi_s));
energyi_v_multiple_run = zeros(length(A(:,1)),length(energyi_s));
vmax_sv_multiple_run = zeros(length(A(:,1)),length(vmax_sv));
vmax_s_multiple_run = zeros(length(A(:,1)),length(vmax_s));
vmax_v_multiple_run = zeros(length(A(:,1)),length(vmax_v));
for vv = 1:length(A(:,1))
[~,~,~,~,...
    vmax_s_multiple_run(vv,:),vmax_v_multiple_run(vv,:),vmax_sv_multiple_run(vv,:),~,~,~,...
    energyi_s_multiple_run(vv,:),energyi_v_multiple_run(vv,:),energyi_sv_multiple_run(vv,:)] = objective_function(...
    nbd,A(vv,1:end-1),TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv);
end

%% solve for mean and std
vmax_s_multiple_run_mean = mean(vmax_s_multiple_run);
vmax_s_multiple_run_std = std(vmax_s_multiple_run);
vmax_v_multiple_run_mean = mean(vmax_v_multiple_run);
vmax_v_multiple_run_std = std(vmax_v_multiple_run);
vmax_sv_multiple_run_mean = mean(vmax_sv_multiple_run);
vmax_sv_multiple_run_std = std(vmax_sv_multiple_run);
energyi_sv_multiple_run_mean = mean(-log(energyi_sv_multiple_run));
energyi_sv_multiple_run_std = std(-log(energyi_sv_multiple_run));
energyi_s_multiple_run_mean = mean(-log(energyi_s_multiple_run));
energyi_s_multiple_run_std = std(-log(energyi_s_multiple_run));
energyi_v_multiple_run_mean = mean(-log(energyi_v_multiple_run));
energyi_v_multiple_run_std = std(-log(energyi_v_multiple_run));

%% plot vmax

rgbColor = hex2rgb(hexColor);
%%%%%%% plot vmax bar graph ps
figure("position",[81,252,1733,512]); % for ter E and Vmax, subplot(1,3,i)
subplot(1,3,2)
X = categorical({'000','001','010','011','100','101','110','111'});
Y = repelem(vmax_s_multiple_run_mean/vmax_s_multiple_run_mean(1),8) ;

b = barh(X,Y);
colororder(rgbColor(2,:))

hold on
er = errorbar(Y, X, repelem(vmax_s_multiple_run_std./vmax_s_multiple_run_mean(1),8), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

xlabel('Vmax')
ylabel('Configuration')
set(gca,'FontSize',20)
title('Ps only (thermo)')
xline(1,'--')
xlim([0 2])

%%%%%%%%%%% plot vmax on bar graph pv
norm_vmax_v_multiple_run_mean = vmax_v_multiple_run_mean./vmax_s_multiple_run_mean(1);
% norm_vmax_v_multiple_run_mean(6) = norm_vmax_v_multiple_run_mean(6)/vmax_v_multiple_run_mean(1);
norm_vmax_v_multiple_run_std = vmax_v_multiple_run_std./vmax_s_multiple_run_mean(1);
% norm_vmax_v_multiple_run_std(6) = norm_vmax_v_multiple_run_std(6)/vmax_v_multiple_run_mean(1);
subplot(1,3,1)
X = categorical({'000','001','010','011','100','101','110','111'});
Y = repelem(norm_vmax_v_multiple_run_mean,8);
b = barh(X,Y);
colororder(rgbColor(2,:))

hold on

er = errorbar(Y, X, repelem(norm_vmax_v_multiple_run_std,8), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

ylabel('Configuration')
xlabel('Vmax')
set(gca,'FontSize',20)
title('Pv only (thermo)')
xline(1,'--')
xlim([0 2])

%%%%%%%%%%% plot vmax on bar graph pspv
norm_vmax_sv_multiple_run_mean = vmax_sv_multiple_run_mean./vmax_s_multiple_run_mean(1);
% norm_vmax_v_multiple_run_mean(6) = norm_vmax_v_multiple_run_mean(6)/vmax_v_multiple_run_mean(1);
norm_vmax_sv_multiple_run_std = vmax_sv_multiple_run_std./vmax_s_multiple_run_mean(1);
% norm_vmax_v_multiple_run_std(6) = norm_vmax_v_multiple_run_std(6)/vmax_v_multiple_run_mean(1);
subplot(1,3,3)
X = categorical({'000','001','010','011','100','101','110','111'});
Y = repelem(norm_vmax_sv_multiple_run_mean,8);
b = barh(X,Y);
colororder(rgbColor(2,:))

hold on

er = errorbar(Y, X, repelem(norm_vmax_sv_multiple_run_std,8), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

ylabel('Configuration')
xlabel('Vmax')
set(gca,'FontSize',20)
title('PsPv (thermo)')
xline(1,'--')
xlim([0 3])

saveas(gcf,'Pure_thermo_quaE_vmax.jpg')
saveas(gcf,'Pure_thermo_quaE_vmax.fig')
%% plot binding energy parameters 

%%%%%%%%%%%%%%%%%% Ps (kinetic)
figure("position",[81,252,1733,512]); % for ter E, subplot(1,3,i)

subplot(1,3,2)
rgbColor = hex2rgb(hexColor);
colororder([rgbColor(8,:)])

X = categorical({'1','2','3','s','12','13','1s','23','2s','3s','123','12s','13s','23s','123s'});
X = reordercats(X,{'1','2','3','s','12','13','1s','23','2s','3s','123','12s','13s','23s','123s'});
Y = energyi_s_multiple_run_mean;
b = barh(X,Y);

hold on

er = errorbar(Y, X, energyi_s_multiple_run_std, '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

title('Ps only (thermo)')
xlabel('Energy')
set(gca,'FontSize',20)
xlim([-20 10])

%%%%%%%%%%%%%%%%%% Pv (thermo)

subplot(1,3,1)
hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);
colororder([rgbColor(8,:)])

X = categorical({'1','2','3','v','12','13','1v','23','2v','3v','123','12v','13v','23v','123v'});
X = reordercats(X,{'1','2','3','v','12','13','1v','23','2v','3v','123','12v','13v','23v','123v'});
Y = energyi_v_multiple_run_mean;
b = barh(X,Y);

hold on

er = errorbar(Y, X, energyi_v_multiple_run_std, '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

title('Pv only (thermo)')
xlabel('Energy')
set(gca,'FontSize',20)
xlim([-20 10])

%%%%%%%%%%%%%%%%%% PsPv

subplot(1,3,3)
hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);
colororder([rgbColor(8,:)])

X = categorical({'sv','1sv','2sv','3sv','12sv','13sv','23sv'});
X = reordercats(X,{'sv','1sv','2sv','3sv','12sv','13sv','23sv'});
Y = energyi_sv_multiple_run_mean([15,21,24,25,28,29,30]);
b = barh(X,Y);

hold on

er = errorbar(Y, X, energyi_sv_multiple_run_std([15,21,24,25,28,29,30]), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

title('PsPv (thermo)')
xlabel('Energy')
set(gca,'FontSize',20)
xlim([-20 10])
saveas(gcf,'Pure_thermo_quaE_etnergy.jpg')
saveas(gcf,'Pure_thermo_quaE_etnergy.fig')