%% generate prediction for constructs with more than one 0A box
close all
clear
clc
n = 4;
fn = sprintf('aug_EMSA_n%1d_nonorm_singlesitecurvefit_inunitG.txt',n); 
A = readmatrix(fn);
pars_range_G = A(:,1:end-1);
N_sample = 1000;
TF_conc_hor = 0:0.2:2;
TF_conc_ver = TF_conc_hor';

%% use each set of optimized pars from augmentation to compute sim data; repeat for all 1000 sets of optimized pars
pred_0A12 = zeros(length(TF_conc_hor),N_sample); % 11 rows for each conetnration, 20 columns for prediction at each set of energy
pred_0A13 = zeros(length(TF_conc_hor),N_sample);
pred_0A23 = zeros(length(TF_conc_hor),N_sample);
pred_0A123 = zeros(length(TF_conc_hor),N_sample);

for le = 1:N_sample
    sim_data_multi_site = predict_multi0A_using_single_siteE(pars_range_G(le,:),TF_conc_hor,n);
    pred_0A12(:,le) = sim_data_multi_site(:,1);
    pred_0A13(:,le) = sim_data_multi_site(:,2);
    pred_0A23(:,le) = sim_data_multi_site(:,3);
    pred_0A123(:,le) = sim_data_multi_site(:,4);
end

pred_0A12_mean = mean(pred_0A12,2);
pred_0A12_std = std(pred_0A12,0,2);
pred_0A13_mean = mean(pred_0A13,2);
pred_0A13_std = std(pred_0A13,0,2);
pred_0A23_mean = mean(pred_0A23,2);
pred_0A23_std = std(pred_0A23,0,2);
pred_0A123_mean = mean(pred_0A123,2);
pred_0A123_std = std(pred_0A123,0,2);

pred_bound = [pred_0A12_mean,pred_0A13_mean,pred_0A23_mean,pred_0A123_mean];
pred_bound_std = [pred_0A12_std,pred_0A13_std,pred_0A23_std,pred_0A123_std];

%% plot model predicted two site graph
load('best_fit_model_prediction_colors.mat')
rgbColor = hex2rgb(hexColor);
load('raw_EMSA_data_Dec_2023.mat')
plot_data =  [A12_present(:,2),A13_present(:,2),A23_present(:,2),A123_present(:,2)]; % for plotting, average
plot_data_std = [A12_present(:,3),A13_present(:,3),A23_present(:,3),A123_present(:,3)];


figure('position',[188,100,1300/2,550])
plot_areaerrorbar_wosim(TF_conc_hor,plot_data,plot_data_std,rgbColor(4:-1:1,:))
plot_all_sim_2(TF_conc_hor,pred_bound,pred_bound_std,rgbColor(4:-1:1,:))
xlabel('[0A~P] (μM)')
ylabel('DNA fraction bound')
set(gca,'FontSize',28)
legend('12','13','23','123','Location','northwest','FontSize',25,'NumColumns',1)
ylim([0 1])
xlim([0 2])
f = gcf;
exportgraphics(f,sprintf('predicted_n%1d_multisite_bindingcurve.tif',n),'Resolution',600)