clear
clc
load('citGspo0A_prediction.mat')
load('Physpankspo0A_prediction.mat')
colors = ['#FF0000';'#ffA000';'#0063dc';'#00BfEB'];
mycolors = hex2rgb(colors);

%%%%%%%% use line width 7.4 for this, and data marker 17
% figure('position',[188,100,1300/2,550])
figure('position',[188,100,1300/2*1.5,550*1.5]); % for augmented manipulation plots only
% plot_areaerrorbar_wosim(2:8,real_data_A./real_data_A(1,:),real_data_A_std./real_data_A(1,:),mycolors(1:2,:)) % normalized
% hold on
% plot_areaerrorbar_wosim(2:8,real_data_H./real_data_H(1,:),real_data_H_std./real_data_H(1,:),mycolors(3:4,:)) % normalized
% hold on
% plot_all_sim(2:8,sim_data_A./sim_data_A(1,:),mycolors(1:2,:)) % normalized
% hold on
% plot_all_sim(2:8,sim_data_H./sim_data_H(1,:),mycolors(3:4,:)) % normalized
plot_areaerrorbar_wosim(2:8,real_data_A,real_data_A_std,mycolors(1:2,:))
hold on
plot_areaerrorbar_wosim(2:8,real_data_H,real_data_H_std,mycolors(3:4,:)) 
hold on
plot_all_sim(2:8,sim_data_A,mycolors(1:2,:)) 
hold on
plot_all_sim(2:8,sim_data_H,mycolors(3:4,:))
ylim([0 6000])
xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
[~, hobj, ~, ~]=legend('','','','','P{\ithy-spank}','P{\itv} none','P{\itcitG}','P{\its} none','Location','Northwest','FontSize',32,'NumColumns',4);
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.9);
set(gca,'FontSize',38)
f = gcf;
exportgraphics(f,'sigmaAH_PcitGPhyspank_Pspo0A_sameR.tif','Resolution',600)
