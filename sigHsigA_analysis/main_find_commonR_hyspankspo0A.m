clear
clc
load('sigA_sigH_real_data.mat') % collected from physpank and citG promtoer
load('new_promoter_activity_single.mat')
spo0A_mean = new_real_data_Pv(1:7,end);
real_data_A = [Physpac_mean,spo0A_mean];
real_data_A_std = [Physpac_st,new_Pv_only_std(1:7,end)];

%% different K, same vmax
fun = @(p) fit_Hill_diffK(p,real_data_A);

nvars = 10;
lb = zeros(1,nvars);
ub = zeros(1,nvars);
lb(1:2) = exp(1); % K1, K2
ub(1:2) = exp(6);
lb(3) = 100; % vmax
ub(3) = 10000;
lb(4:10) = 0; % R
ub(4:10) = 2;

filename = 'spo0A_hypspank_expression_fit_together_diffK.txt';

for pp = 1:20
    [x,fval] = particleswarm(fun,nvars,lb,ub);
    create_parameter_file(filename,x,fval,'a+');
end

% prediction
close all
colors = ['#FF0000';'#ffA000'];
mycolors = hex2rgb(colors);
A = readmatrix(filename);
[~,I] = min(A(:,end));
pars = A(I,1:end-1);
[~,sim_data_A] = fit_Hill_diffK(pars,real_data_A);
save('Physpankspo0A_prediction.mat','real_data_A','sim_data_A','real_data_A_std')
%
figure('position',[188,100,1300/2,550])
plot_areaerrorbar_wosim(2:8,real_data_A,real_data_A_std,mycolors)
plot_all_sim(2:8,sim_data_A,mycolors)
ylim([0 5600])
xlabel('Time (h)')
ylabel('Promoter Activity')
legend('Phy-spank','spo0A Pv','Location','Northwest')
set(gca,'FontSize',22)
f = gcf;
exportgraphics(f,'sigmaH_Physpank_Pspo0A_sameR_diffK.tif','Resolution',600)

%% different vmax, same K
fun = @(p) fit_Hill_diffvmax(p,real_data_A);

nvars = 10;
lb = zeros(1,nvars);
ub = zeros(1,nvars);
lb(1) = exp(1); % K
ub(1) = exp(6);
lb(2:3) = 100; % vmax1, vmax2
ub(2:3) = 10000;
lb(4:10) = 0; % R
ub(4:10) = 2;

filename = 'spo0A_hypspank_expression_fit_together_diffvmax.txt';

for pp = 1:20
    [x,fval] = particleswarm(fun,nvars,lb,ub);
    create_parameter_file(filename,x,fval,'a+');
end

% prediction
close all
colors = ['#FF0000';'#ffA000'];
mycolors = hex2rgb(colors);
A = readmatrix(filename);
[~,I] = min(A(:,end));
pars = A(I,1:end-1);
[~,sim_data_A] = fit_Hill_diffvmax(pars,real_data_A);
save('Physpankspo0A_prediction.mat','real_data_A','sim_data_A','real_data_A_std')
%
figure('position',[188,100,1300/2,550])
plot_areaerrorbar_wosim(2:8,real_data_A,real_data_A_std,mycolors)
plot_all_sim(2:8,sim_data_A,mycolors)
ylim([0 5600])
xlabel('Time (h)')
ylabel('Promoter Activity')
legend('Phy-spank','spo0A Pv','Location','Northwest')
set(gca,'FontSize',22)
f = gcf;
exportgraphics(f,'sigmaH_Physpank_Pspo0A_sameR_diffvmax.tif','Resolution',600)