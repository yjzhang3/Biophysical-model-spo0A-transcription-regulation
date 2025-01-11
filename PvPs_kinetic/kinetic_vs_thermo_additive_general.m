%% replicate paper figure under the parameter randomly sampled
clear
load('thermo_general_not_additive_example_pars.mat')
load('kinetic_always_additive_example_pars.mat')

nbd = 5; % two promoters and three 0A boxes
mut_mat = [0,1,1]; % there is only one 0A box unmutated on DNA fragment
mut_mat_sv = mut_mat;
TF_conc_t = (get_aps(2:9));
H_conc_t = zeros(1,length(TF_conc_t))+0.01;
A_conc_t = zeros(1,length(TF_conc_t))+0.01; % holoenzyme concentration same as used for sampling
group_array_s.g1 = 1:8;
group_array_v.g1 = 1:8;
group_array_sv.g1 = 1:8; % always true for thermo model
vmax_array_s.g1 = 5;
vmax_array_v.g1 = 5; % same as what it's used for samplign
vmax_array_sv.g1 = vmax_array_s.g1 + vmax_array_v.g1;

sim_data_Ps_rep_thermo = time_dep_TR_new_wSigma( ...
        nbd-1,energyi_s_thermo,mut_mat,TF_conc_t, ...
        H_conc_t,group_array_s,vmax_array_s);
    
sim_data_Ps_rep_kinetic = time_dep_TR_new_wSigma( ...
        nbd-1,energyi_s_kinetic,mut_mat,TF_conc_t, ...
        H_conc_t,group_array_s,vmax_array_s);

sim_data_Pv_rep_thermo = time_dep_TR_new_wSigma( ...
        nbd-1,energyi_v_thermo,mut_mat,TF_conc_t, ...
        A_conc_t,group_array_v,vmax_array_v);

sim_data_Pv_rep_kinetic = time_dep_TR_new_wSigma( ...
        nbd-1,energyi_v_kinetic,mut_mat,TF_conc_t, ...
        A_conc_t,group_array_v,vmax_array_v);

sim_data_PsPv_rep_thermo = time_dep_TR_new_wSigma_twoP( ...
    nbd,energyi_sv_thermo,mut_mat_sv,TF_conc_t,...
    H_conc_t,A_conc_t,...
    vmax_array_s,group_array_s,...
    vmax_array_v,group_array_v,...
    vmax_array_sv,group_array_sv);

sim_data_PsPv_rep_kinetic = time_dep_TR_new_wSigma_twoP( ...
    nbd,energyi_sv_kinetic,mut_mat_sv,TF_conc_t,...
    H_conc_t,A_conc_t,...
    vmax_array_s,group_array_s,...
    vmax_array_v,group_array_v,...
    vmax_array_sv,group_array_sv);

%% load real data
mut_mat_all = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
strain_type = find(ismember(mut_mat_all,mut_mat,'row')); % choose the right color to plot
hexColor=['#000000';'#4CEB4E';'#e6d800';'#9b19f5';'#39F9D9';'#ffa300';'#dc0ab4';'#0bb4ff'];
rgbColor_raw = hex2rgb(hexColor);
rgbColor = rgbColor_raw(strain_type,:);

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

ratio_real = (real_data_Ps+real_data_Pv)./real_data_PsPv;
ratio_std = (real_data_Ps_std+real_data_Pv_std)./real_data_PsPv;

% lgd_add = {'','WT','','1*','','2*','','3*','','12*','','13*','','23*','','123*'};

%% plot
clear lgd
figure('position',[188,100,1300/2,550])
plot_areaerrorbar_wosim(ratio_real,ratio_std,rgbColor)
ratio_thermo = (sim_data_Pv_rep_thermo+sim_data_Ps_rep_thermo)./sim_data_PsPv_rep_thermo;
ratio_kinetic = (sim_data_Pv_rep_kinetic+sim_data_Ps_rep_kinetic)./sim_data_PsPv_rep_kinetic;
hold on
plot(2:9,ratio_thermo,'LineWidth',3.78,'Color',rgbColor,'LineStyle','--','Marker','x','MarkerSize',20)
hold on
plot(2:9,ratio_kinetic,'LineWidth',3.78,'Color',rgbColor,'LineStyle','--','Marker','o','MarkerSize',15)
    
ylim([0 2])
xlim([2 9])
lgd_add_sim_t = {'Simulated 123 (thermodynamic)','Simulated 23 (thermodynamic)','Simulated 13 (thermodynamic)','Simulated 12 (thermodynamic)','Simulated 3 (thermodynamic)','Simulated 2 (thermodynamic)','Simulated 1 (thermodynamic)','Simulated none (thermodynamic)'};
lgd_add_sim_k = {'Simulated 123 (kinetic)','Simulated 23 (kinetic)','Simulated 13 (kinetic)','Simulated 12 (kinetic)','Simulated 3 (kinetic)','Simulated 2 (kinetic)','Simulated 1 (kinetic)','Simulated none (kinetic)'};
lgd_add_real = {'','Experimental WT','','Experimental 23','','Experimental 13','','Experimental 12','','Experimental 3','','Experimental 2','','Experimental 1','','Experimental none'};
A = strain_type*2-1;
B = A+1;
AB = reshape([A;B], size(A,1), []);
lgd{1} = [lgd_add_real(AB(1))];
lgd{2} = [lgd_add_real(AB(2))];
lgd{3} = lgd_add_sim_t(strain_type);
lgd{4} = lgd_add_sim_k(strain_type);
legend(string(lgd),'NumColumns',1,'Location','northeast','FontSize',22)
xlabel('Time (h)')
ylabel('Additivity Level')
set(gca,'FontSize',28)
% title('Purely Thermodynamic')

% saveas(gcf,'kinetic_vs_thermo_additive_general.fig')
f = gcf;
exportgraphics(f,'kinetic_vs_thermo_additive_general.tif','Resolution',300)