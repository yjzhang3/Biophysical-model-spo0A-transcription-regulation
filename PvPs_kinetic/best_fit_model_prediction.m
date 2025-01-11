%% this file reads best parameters from optimization and plot vmax and energy 
clc
clear
close all

%% 
fn = 'purekinetic_Sept_2024_extrarun_highera.txt';
A = readmatrix(fn);
[M,I] = min(A(:,end));
pars = A(I,:);
pars = pars(1:end-1);

%% concentration  
TF_conc_t = get_aps(2:9);
TF_conc_const = zeros(1,8)+TF_conc_t(1);
nbd = 5;
nbd_s = 4;
nbd_v = 4;
mut_mat_s = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_v = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
mut_mat_sv = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];

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

%% vmax grouping
group_array_s.g1 = 1; 
group_array_s.g2 = 2; 
group_array_s.g3 = 3; 
group_array_s.g4 = 4; 
group_array_s.g5 = 5; 
group_array_s.g6 = 6; 
group_array_s.g7 = 7; 
group_array_s.g8 = 8; 

group_array_v.g1 = 1; 
group_array_v.g2 = 2; 
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

%% generate simulated data
[final_diff,sim_data_s,sim_data_v,sim_data_sv,vmax_s_final,vmax_v_final,vmax_sv_final,vmax_array_s,vmax_array_v,vmax_array_sv,...
    energyi_s,energyi_v,energyi_sv,H_conc_t,A_conc_t] = objective_function(...
    nbd,pars,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv);
save('sim_data_original.mat','sim_data_s','sim_data_v','sim_data_sv')
save('best_fit_energy.mat','energyi_s','energyi_sv','energyi_v')
save('best_fit_vmax.mat','vmax_s_final','vmax_v_final','vmax_array_s','vmax_array_v','vmax_array_sv','group_array_s','group_array_v')
save('best_fit_RNAP_t.mat','H_conc_t','A_conc_t')

%%
% sim data with constant RNA pol
% [~,sim_data_s_constsig,sim_data_v_constsig,sim_data_sv_constsig,~,~,~,~,~,~,...
%     ~,~,~,H_conc_t_norm,A_conc_t_norm] = objective_function(...
%     nbd,pars,TF_conc_t,...
%     mut_mat_s,mut_mat_v,mut_mat_sv, ...
%     real_data_Ps,real_data_Pv,real_data_PsPv,...
%     group_array_s,group_array_v,group_array_sv);
% save('sim_data_const_holoenzyme.mat','sim_data_s_constsig','sim_data_v_constsig','sim_data_sv_constsig')
%% 
% sim data with constant TF but nonconstant RNA pol
% [diff,sim_data_s_constTF,sim_data_v_constTF,sim_data_sv_constTF,~,~,~,~,~,~,...
%     ~,~,~,H_conc_t,A_conc_t] = objective_function(...
%     nbd,pars,TF_conc_const,...
%     mut_mat_s,mut_mat_v,mut_mat_sv, ...
%     real_data_Ps,real_data_Pv,real_data_PsPv,...
%     group_array_s,group_array_v,group_array_sv);
% save('sim_data_const_0AP.mat','sim_data_s_constTF','sim_data_v_constTF','sim_data_sv_constTF')


%% obtain results for multiple run
% rows = find(A(:,end)<0.04);
% A_final = A(rows,:);

A_final = A;
energyi_sv_multiple_run = zeros(length(A_final(:,1)),length(energyi_sv));
energyi_s_multiple_run = zeros(length(A_final(:,1)),length(energyi_s));
energyi_v_multiple_run = zeros(length(A_final(:,1)),length(energyi_v));
vmax_s_multiple_run = zeros(length(A_final(:,1)),length(vmax_s_final));
vmax_v_multiple_run = zeros(length(A_final(:,1)),length(vmax_v_final));
vmax_sv_multiple_run = zeros(length(A_final(:,1)),length(vmax_sv_final));
H_conc_t_multiple_run = zeros(length(A_final(:,1)),length(H_conc_t));
A_conc_t_multiple_run = zeros(length(A_final(:,1)),length(A_conc_t));
for vv = 1:length(A_final(:,1))
    [~,~,~,~,...
        vmax_s_multiple_run(vv,:),vmax_v_multiple_run(vv,:),vmax_sv_multiple_run(vv,:),~,~,~,...
        energyi_s_multiple_run(vv,:),energyi_v_multiple_run(vv,:),energyi_sv_multiple_run(vv,:),...
        H_conc_t_multiple_run(vv,:),A_conc_t_multiple_run(vv,:)] = ...
    objective_function(...
    nbd,A_final(vv,1:end-1),TF_conc_t,...
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
H_conc_t_multiple_run_mean = mean(H_conc_t_multiple_run);
H_conc_t_multiple_run_std = std(H_conc_t_multiple_run);
A_conc_t_multiple_run_mean = mean(A_conc_t_multiple_run);
A_conc_t_multiple_run_std = std(A_conc_t_multiple_run);

%% plot sigma dynamics predicted
colors = ['#ffA000';'#00BfEB'];
mycolors = hex2rgb(colors);
figure('position',[188,100,1300/2,550])
plot_all_sim_2(2:9,H_conc_t_multiple_run_mean',H_conc_t_multiple_run_std',mycolors(1,:))
hold on
plot_all_sim_2(2:9,A_conc_t_multiple_run_mean',A_conc_t_multiple_run_std',mycolors(2,:))
xlabel('Time (h)')
ylabel('Concentration (μM)')
set(gca,'FontSize',28)
xlim([2 9])
ylim([0 2.5])
legend('','RNAP-\sigma^{H}','','RNAP-\sigma^{A}','location','northwest','FontSize',22)
f = gcf;


f = gcf;
exportgraphics(f,'SigmaDynamics_WT.tif','Resolution',600)

%% plot vmax
% 
load('best_fit_model_prediction_colors.mat')

%%%%%%% plot vmax bar graph ps
figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
X = categorical({'000','100','010','001','110','101','011','111'});
X = reordercats(X,{'000','100','010','001','110','101','011','111'});
Y_s = [vmax_s_multiple_run_mean(8),vmax_s_multiple_run_mean(1:7)];
Y_s = Y_s./vmax_s_multiple_run_mean(8);


% subplot(1,3,3)
b = barh(X,Y_s);
set(gca,'Ydir','reverse')
colororder([0.5 0.5 0.5])

hold on
er = errorbar(Y_s, X,vmax_s_multiple_run_std([8,1:7])./vmax_s_multiple_run_mean(8), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

xlabel('v_{max}')
ylabel('Configuration')
set(gca,'FontSize',32)
% title('Ps only (kinetic)')
xline(1,'--')
xlim([0 2])
f = gcf;
exportgraphics(f,'Pure_kinetic_Ps_vmax.tif','Resolution',600)

%%%%%%%%%%% plot vmax on bar graph pv

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
X = categorical({'000','100','010','001','110','101','011','111'});
X = reordercats(X,{'000','100','010','001','110','101','011','111'});
% Y = [1,vmax_v(1:4),vmax_v(1)*vmax_v(3),vmax_v(5:end)];
Y_v = [vmax_v_multiple_run_mean(8),vmax_v_multiple_run_mean(1:7)];
Y_v = Y_v./vmax_v_multiple_run_mean(8);
b = barh(X,Y_v);
set(gca,'Ydir','reverse')
colororder([0.5 0.5 0.5])
hold on

er = errorbar(Y_v, X, vmax_v_multiple_run_std([8,1:7])./vmax_v_multiple_run_mean(8), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;
ylabel('Configuration')
xlabel('v_{max}')
set(gca,'FontSize',32)
% title('Pv only (kinetic)')
xline(1,'--')
xlim([0 2])
f = gcf;
exportgraphics(f,'Pure_kinetic_Pv_vmax.tif','Resolution',600)

%%%%%%%%%%% plot vmax on bar graph pvps

figure('position',[188,100,1300/2,550]); % for strain fits, single row subplot(2,1,1)
X = categorical({'000','100','010','001','110','101','011','111'});
X = reordercats(X,{'000','100','010','001','110','101','011','111'});
% Y = [1,vmax_v(1:4),vmax_v(1)*vmax_v(3),vmax_v(5:end)];
Y = vmax_sv_multiple_run_mean([8,1:7])./vmax_sv_multiple_run_mean(8);
b = barh(X,Y);
set(gca,'Ydir','reverse')
colororder([0.5 0.5 0.5])
hold on

er = errorbar(Y, X, vmax_sv_multiple_run_std([8,1:7])./vmax_sv_multiple_run_mean(8), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;
ylabel('Configuration')
xlabel('v_{max}')
set(gca,'FontSize',32)
% title('PvPs (kinetic)')
xline(1,'--')
xlim([0 2])

f = gcf;
exportgraphics(f,'Pure_kinetic_PvPs_vmax.tif','Resolution',600)
%% plot binding energy parameters 

%%%%%%%%%%%%%%%%%% Ps (kinetic)
figure("position",[81,252,1733,512]); % for ter E, subplot(1,3,i)

subplot(1,3,3)
% hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
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

title('Ps only (kinetic)')
xlabel('Energy')
set(gca,'FontSize',20)
xlim([-20 10])

%%%%%%%%%%%%%%%%%% Pv (thermo)

subplot(1,3,2)
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

subplot(1,3,1)
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

title('PsPv (kinetic)')
xlabel('Energy')
set(gca,'FontSize',20)
xlim([-20 10])

% saveas(gcf,'Pure_kinetic_energy.fig')
saveas(gcf,'Pure_kinetic_energy.tif')

%% since interaction energies and binding energies are shared by two promoters, just plot these energies
%%%%%%%%%%%%%%%%%% Ps (kinetic)
figure("position",[81,252,1733/3,512]); % for ter E, subplot(1,3,i)

% subplot(1,3,3)
% hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);
colororder([rgbColor(8,:)])

X = categorical({'G_{12}','G_{13}','G_{23}','G_{123}'});
X = reordercats(X,{'G_{12}','G_{13}','G_{23}','G_{123}'});
Y = energyi_s_multiple_run_mean([5,6,8,11]);
b = barh(X,Y);

hold on

er = errorbar(Y, X, energyi_s_multiple_run_std([5,6,8,11]), '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

% title('Ps only (kinetic)')
xlabel('Energy')
set(gca,'FontSize',20)
xlim([-20 10])
