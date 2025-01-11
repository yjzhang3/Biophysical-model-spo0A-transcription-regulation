clear
close all
load('sim_data_const_holoenzyme.mat')

%%
figure("position",[39,179,1800,636]);
subplot(1,3,1)
load('best_fit_model_prediction_colors.mat')

mycolors = hex2rgb(hexColor);
strain_name = {'123','23','13','12','3','2','1','none'};
sim_data = sim_data_sv_constsig;
for i = [7:-1:2,1,8]
    plot(2:9,sim_data(:,i),'LineWidth',2.78,'LineStyle','--','Marker',"diamond",'MarkerFaceColor',mycolors(i,:),'Color',mycolors(i,:))
    hold on
    % predicted 123* normalized by RNAp
%     plot(2:9,sim_data(:,8),'LineWidth',2,'LineStyle','--','Marker',"diamond",'Color',mycolors(8,:))
    xlabel('Time (h)')
    ylabel('v_{MS}^{effective} (a.u.)')
    set(gca,'FontSize',22)
    xlim([2 9])
    ylim([0 200])
end
% title('PsPv strains')
legend(string(strain_name([7:-1:2,1,8])),'FontSize',20,'NumColumns',2,'Location','southwest')


% figure();
subplot(1,3,2)
mycolors = hex2rgb(hexColor);
sim_data = sim_data_v_constsig;
for i = [7:-1:2,1,8]
%     subplot(7,1,i)
    plot(2:9,sim_data(:,i),'LineWidth',2.78,'LineStyle','--','Marker',"diamond",'MarkerFaceColor',mycolors(i,:),'Color',mycolors(i,:))
    hold on
    % predicted 123* normalized by RNAp
%     plot(2:9,sim_data(:,8),'LineWidth',2,'LineStyle','--','Marker',"diamond",'Color',mycolors(8,:))
    xlabel('Time (h)')
    ylabel('v_{MS}^{effective} (a.u.)')
 set(gca,'FontSize',22)
    xlim([2 9])
    ylim([0 200])
%     hold on
end
% title('Pv only strains')
legend(string(strain_name([7:-1:2,1,8])),'FontSize',20,'NumColumns',2,'Location','southwest')
% f = gcf;
% exportgraphics(f,'Constant_sigma_allstrains_sim_Pv.tif','Resolution',600)
% fc = (sim_data(end,:)-sim_data(end,end))./sim_data(end,end);

% figure();
subplot(1,3,3)
sim_data = sim_data_s_constsig;
for i = [7:-1:2,1,8]
%     subplot(7,1,i)
    plot(2:9,sim_data(:,i),'LineWidth',2.78,'LineStyle','--','Marker',"diamond",'MarkerFaceColor',mycolors(i,:),'Color',mycolors(i,:))
    hold on
    % predicted 123* normalized by RNAp
%     plot(2:9,sim_data(:,8),'LineWidth',2,'LineStyle','--','Marker',"diamond",'Color',mycolors(8,:))
    xlabel('Time (h)')
    ylabel('v_{MS}^{effective} (a.u.)')
    set(gca,'FontSize',22)
    xlim([2 9])
    ylim([0 200])
%     hold on
end
% title('Ps only strains')
fc = (sim_data(end,:)-sim_data(end,end))./sim_data(end,end);

legend(string(strain_name([7:-1:2,1,8])),'FontSize',20,'NumColumns',2,'Location','northwest')
f = gcf;
exportgraphics(f,'Constant_sigma_allstrains_sim.tif','Resolution',600)
