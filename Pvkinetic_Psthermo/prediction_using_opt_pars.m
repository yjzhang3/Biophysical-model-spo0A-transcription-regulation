clear
clc
close all
% fn = 'July7-2024-Pvkinetic-Psthermo-quaE-noplugin-wEMSA-nosv-highera.txt';
fn = 'pvkinetic_psthermo_Sept_2024_extrarun_highera.txt';
A = readmatrix(fn);
[M,I] = min(A(:,end));
pars = A(I,:);
pars = pars(1:end-1);
% pars = [0.213247,1.000000,1.000000,22026.434615,0.000045,2.131833,1.000000,... % beta, s, v, 12, 13, 1s,1v
%     22026.465795,0.154449,1.000000,... % 23, 2s, 2v
%     0.659886,1.000000,1.000000,... % 3s, 3v, sv
%     1,1,1,1,1,1,...   % 123,12s,12v,13s,13v,1sv
%     1,1,1,1,... % 23s, 23v, 2sv, 3sv
%     1,1,1,1,1,...% 123s, 123v, 12sv, 13sv, 23sv
%     1181.880436,348.145215,693.033560,497.376537,315.533014,578.809233,294.407042,10.027711,... % Vmaxs000, vmaxv1-7
%     0.052233,0.123270,0.515185,0.934759,1.208462,1.517096,1.999573,1.999998];
% pars = [0.328983,1.000000,1.000000,148.413159,0.006738,1.648092,1.000000,148.413159,0.100392,1.000000,0.709883,1.000000,1.000000,0.584680,1.253634,1.000000,1.101513,1.000000,1.000000,0.801655,1.000000,1.000000,1.000000,0.951701,1.000000,1.000000,1.000000,1.000000,1342.687267,443.260979,729.757633,605.294161,276.892759,519.710142,269.441680,10.030857,679.911123,0.046564,0.149758,0.401169,0.731033,1.581299,1.157686,1.999637,1.999998,0.334836,0.383666,0.461009,0.584288,1.104898,0.977449,1.660719,1.840029,0.357165];

%%
group_array_s.g1 = 1:8;
group_array_v.g1 = 1; %000
group_array_v.g2 = 2; %001
group_array_v.g3 = 3; 
group_array_v.g4 = 4; 
group_array_v.g5 = 5; 
group_array_v.g6 = 6; 
group_array_v.g7 = 7; 
group_array_v.g8 = 8; 


group_array_sv.g1 = 1; 
group_array_sv.g2 = 2; 
group_array_sv.g3 = 3; 
group_array_sv.g4 = 4; 
group_array_sv.g5 = 5; 
group_array_sv.g6 = 6; 
group_array_sv.g7 = 7; 
group_array_sv.g8 = 8; 


load('new_promoter_activity_single.mat');
load('new_PsPv_promoter_activity.mat');

real_data_Ps = new_real_data_Ps(1:end-1,:); % exclude the real WT, instead WT is Ps only, no mutations
% real_data_Ps = Ps_rand;
real_data_Ps_std = new_Ps_only_std(1:end-1,:);

real_data_Pv = new_real_data_Pv(1:end-1,:);
% real_data_Pv = Pv_rand;
real_data_Pv_std = new_Pv_only_std(1:end-1,:);

real_data_PsPv = PvPs_avg_new(1:end-1,:);

TF_conc_t = get_aps(2:9);

% nbd = 6;
nbd = 5;
% mut_mat_s = [[1,1,1,1];[0,1,1,1];[1,0,1,1];[1,1,0,1];[0,0,1,1];[0,1,0,1];[1,0,0,1];[0,0,0,1]];
mut_mat_s = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_v = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
% mut_mat_sv = [[1,1,1,1];[0,1,1,1];[1,0,1,1];[1,1,0,1];[0,0,1,1];[0,1,0,1];[1,0,0,1];[0,0,0,1]];
mut_mat_sv = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];

%% generate sim data

% energy related pars
% beta,s,v, (1-3)
% 12,13,1s,1v, (4-7)
% 23,2s,2v, (8-10)
% 3s,3v, (11-12)
% sv (13)

% pars([4:5,8]) = [-3,-0.3,-3];
% pars([6,9,11]) = [40,30,-3];

[final_diff,sim_data_Ps,sim_data_Pv,sim_data_PsPv,vmax_s_final,vmax_v_final,vmax_sv_final,H_conc_t,A_conc_t,...
    energyi_s,energyi_v,energyi_sv] = objective_function(...
    nbd,pars,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv);

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
load('best_fit_model_prediction_colors.mat')
rgbColor = hex2rgb(hexColor);

% unnorm,type 1234
figure('position',[188,174,1177,510]); % for strain fits, single row
% subplot(2,2,1)
subplot(1,2,1)
plot_areaerrorbar_wsim(real_data_Ps(:,1:4),real_data_Ps_std(:,1:4),sim_data_Ps(:,1:4),rgbColor(1:4,:));
title('Ps-only, Single Mutations')
legend('','Real WT','Sim WT','','Real 1*','Sim 1*','','Real 2*','Sim 2*','','Real 3*','Sim 3*','Location','northwest','NumColumns',2,'FontSize',16)
xlabel('Time (hours)')
ylabel('Transcription Rate')
ylim([0 1200])
xlim([2 9])
set(gca,'FontSize',20)

% unnorm,type 15678
% subplot(2,2,2)
subplot(1,2,2)
plot_areaerrorbar_wsim(real_data_Ps(:,[1,5:8]),real_data_Ps_std(:,[1,5:8]),sim_data_Ps(:,[1,5:8]),rgbColor([1,5:8],:))
title('Ps-only, Multiple Mutations')
legend('','Real WT','Sim WT','','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*','Location','northwest','NumColumns',2,'FontSize',12)
xlabel('Time (hours)')
ylabel('Transcription Rate')
ylim([0 1200])
xlim([2 9])
set(gca,'FontSize',20)
saveas(gcf,'PvkineticPsthermo_Ps.fig')
saveas(gcf,'PvkineticPsthermo_Ps.jpg')

%% Pv real and simulation on the same graph unnormalized and normalized
figure('position',[188,174,1177,510]); % for strain fits, single row

% subplot(2,2,1)
subplot(1,2,1)
plot_areaerrorbar_wsim(real_data_Pv(:,1:4),real_data_Pv_std(:,1:4),sim_data_Pv(:,1:4),rgbColor(1:4,:));
title('Pv-only, Single Mutations')
legend('','Real WT','Sim WT','','Real 1*','Sim 1*','','Real 2*','Sim 2*','','Real 3*','Sim 3*','Location','northwest','NumColumns',2,'FontSize',16)
xlabel('Time (hours)')
ylabel('Transcription Rate')
ylim([0 500])
xlim([2 9])
set(gca,'FontSize',20)

% subplot(2,2,2)
subplot(1,2,2)
plot_areaerrorbar_wsim(real_data_Pv(:,[1,5:8]),real_data_Pv_std(:,[1,5:8]),sim_data_Pv(:,[1,5:8]),rgbColor([1,5:8],:))
title('Pv-only, Multiple Mutations')
legend('','Real WT','Sim WT','','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*','NumColumns',2,'Location','northwest','FontSize',12)
xlabel('Time (hours)')
ylabel('Transcription Rate')
ylim([0 500])
xlim([2 9])
set(gca,'FontSize',20)
saveas(gcf,'PvkineticPsthermo_Pv.fig')
saveas(gcf,'PvkineticPsthermo_Pv.jpg')
%% see two promoter model fits to pspv data well
rgbColor = hex2rgb(hexColor);
load('new_PsPv_promoter_activity.mat')

sim_all = sim_data_PsPv;

figure('position',[188,174,1177,510]); % for strain fits, single row
subplot(1,2,1)
plot_areaerrorbar_wsim(real_data_PsPv(:,1:4),PvPs_std_new(1:end-1,1:4),sim_all(:,1:4),rgbColor(1:4,:))
legend('','Real WT','Sim WT','','Real 1*','Sim 1*','','Real 2*','Sim 2*','','Real 3*','Sim 3*','Location','northwest','NumColumns',2,'FontSize',16)
title('PsPv Together, Single Mutation')
xlabel('time')
ylabel('Transcription Rate')
xlim([2 9])
ylim([0 1400])
set(gca,'FontSize',20)

subplot(1,2,2)
plot_areaerrorbar_wsim(PvPs_avg_new(1:end-1,[1,5:8]),PvPs_std_new(1:end-1,[1,5:8]),sim_all(:,[1,5:8]),rgbColor([1,5:8],:))
title('PsPv Together, Multiple Mutations')
legend('','Real WT','Sim WT','','Real 12*','Sim 12*','','Real 13*','Sim 13*','','Real 23*','Sim 23*','','Real 123*','Sim 123*','Location','northwest','NumColumns',2,'FontSize',12)
xlabel('time')
ylabel('Transcription Rate')
xlim([2 9])
ylim([0 1400])
set(gca,'FontSize',20)
% sgtitle('PvPs-together')

saveas(gcf,'PvkineticPsthermo_PsPv.fig')
saveas(gcf,'PvkineticPsthermo_PsPv.jpg')

%% obtain vmax and energy for multiple runs
energyi_sv_multiple_run = zeros(length(A(:,1)),length(energyi_sv));
energyi_s_multiple_run = zeros(length(A(:,1)),length(energyi_s));
energyi_v_multiple_run = zeros(length(A(:,1)),length(energyi_v));
vmax_sv_multiple_run = zeros(length(A(:,1)),length(vmax_sv_final));
vmax_s_multiple_run = zeros(length(A(:,1)),length(vmax_s_final));
vmax_v_multiple_run = zeros(length(A(:,1)),length(vmax_v_final));
for vv = 1:length(A(:,1))
[~,~,~,~,...
    vmax_s_multiple_run(vv,:),vmax_v_multiple_run(vv,:),vmax_sv_multiple_run(vv,:),~,~,...
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

hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

%%%%%%% plot vmax bar graph ps
figure("position",[81,252,1733,512]); % for ter E and Vmax, subplot(1,3,i)
X = categorical({'000','001','010','011','100','101','110','111'});
Y_s = repelem(vmax_s_multiple_run_mean(1),8);
Y_s = Y_s./vmax_s_multiple_run_mean(1);

subplot(1,3,2)
b = barh(X,Y_s);
colororder(rgbColor(2,:))

hold on
er = errorbar(Y_s, X, repelem(vmax_s_multiple_run_std./vmax_s_multiple_run_mean(1),8), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

xlabel('Vmax')
ylabel('Configuration')
set(gca,'FontSize',20)
title('Ps only (thermo)')
xline(1,'--')
xlim([0 2])

%%%%%%%%%%% plot vmax on bar graph pv

subplot(1,3,1)
X = categorical({'000','001','010','011','100','101','110','111'});
% X = reordercats(X,{'vmax_{000}','vmax_{100}','vmax_{010}','vmax_{001}','vmax_{110}','vmax_{011}','vmax_{111}'});
% Y = [1,vmax_v(1:4),vmax_v(1)*vmax_v(3),vmax_v(5:end)];
Y_v = vmax_v_multiple_run_mean; 
Y_v = Y_v./vmax_s_multiple_run_mean(1);
b = barh(X,Y_v);
colororder(rgbColor(2,:))

hold on

er = errorbar(Y_v, X, vmax_v_multiple_run_std./vmax_s_multiple_run_mean(1), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;
ylabel('Configuration')
xlabel('Vmax')
set(gca,'FontSize',20)
title('Pv only (kinetic)')
xline(1,'--')
xlim([0 2])

%%%%%%%%%%% plot vmax on bar graph pvps

subplot(1,3,3)
X = categorical({'000','001','010','011','100','101','110','111'});
% X = reordercats(X,{'vmax_{000}','vmax_{100}','vmax_{010}','vmax_{001}','vmax_{110}','vmax_{011}','vmax_{111}'});
% Y = [1,vmax_v(1:4),vmax_v(1)*vmax_v(3),vmax_v(5:end)];
Y = vmax_sv_multiple_run_mean./vmax_s_multiple_run_mean(1);
b = barh(X,Y);
colororder(rgbColor(2,:))

hold on

er = errorbar(Y, X, vmax_sv_multiple_run_std./vmax_s_multiple_run_mean(1), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;
ylabel('Configuration')
xlabel('Vmax')
set(gca,'FontSize',20)
title('PsPv (mixed)')
xline(1,'--')
xlim([0 3])

saveas(gcf,'PvkineticPsthermo_vmax.fig')
saveas(gcf,'PvkineticPsthermo_vmax.jpg')

%% plot binding energy parameters 

%%%%%%%%%%%%%%%%%% Ps (kinetic)
figure("position",[81,252,1733,512]); % for ter E, subplot(1,3,i)

subplot(1,3,2)
hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
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

title('Pv only (kinetic)')
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

title('PsPv (mixed)')
xlabel('Energy')
set(gca,'FontSize',20)
xlim([-20 10])

saveas(gcf,'PvkineticPsthermo_energy.fig')
saveas(gcf,'PvkineticPsthermo_energy.jpg')

