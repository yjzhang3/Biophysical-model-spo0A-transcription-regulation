function fit_original_dataset(n,range1)
n_exp_per_strain = [4,3,3,3,3,2,2]; % number of experiments collected for each strain

%% data loading and processing

load('raw_EMSA_Data_detail.mat')
load('raw_EMSA_Data_Dec_2023.mat')

% real data format: each n column is the fraction bound for each type of
% construct (for fitting)
real_data = [A1_present_raw,A2_present_raw,A3_present_raw];

% for plotting
plot_data =  [A1_present(:,2),A2_present(:,2),A3_present(:,2)];
plot_data_std = [A1_present(:,3),A2_present(:,3),A3_present(:,3)];

TF_conc = 0:0.2:2;
TF_conc = TF_conc';

fun1 = @(G) objective_G_together_strain_weight(n,G,TF_conc,real_data,n_exp_per_strain);
nvars = 3;
lb = zeros(1,3)+range1(1);
ub = zeros(1,3)+range1(2);

x_all = zeros(20,nvars);
err_all = zeros(20,1);
fn = sprintf('original_EMSA_n%1d_nonorm_singlesitecurvefit.txt',n);
for i = 1:20
[x,err] = particleswarm(fun1,nvars,lb,ub);
x_all(i,:) = x;
err_all(i) = err;
create_parameter_file(fn, x, err, 'a+')
end

%% plot best result when fitted 
load('best_fit_model_prediction_colors.mat')
rgbColor = hex2rgb(hexColor);
[M,I] = min(err_all);
pars = x_all(I,:);

[~,~,sim_plot_data] =  objective_G_together_strain_weight(n,pars,TF_conc,real_data,n_exp_per_strain);

figure('position',[188,100,1300/2,550])
plot_areaerrorbar_wosim(TF_conc,plot_data(:,1:3),plot_data_std(:,1:3),rgbColor(7:-1:5,:))
hold on
plot_all_sim(TF_conc,sim_plot_data(:,1:3),rgbColor(7:-1:5,:));
xlabel('[0A~P] (μM)')
ylabel('DNA fraction bound')
set(gca,'FontSize',28)
% legend('','Real 234*','Sim 234*','','Real 134*','Sim 134*','','Real 124*','Sim 124*','Location','northwest','FontSize',12)
% legend('','None','','','12*','','','13*','','','23*','','Location','northwest','FontSize',32)
legend('1','2','3','Location','northwest','FontSize',25,'NumColumns',1)

ylim([0 1])
xlim([0 2])
f = gcf;
fn = sprintf('original_EMSA_n%1d_nonorm_singlesitecurvefit.tif',n);
exportgraphics(f,fn,'Resolution',600)
