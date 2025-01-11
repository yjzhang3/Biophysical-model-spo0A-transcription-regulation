%% plot tf dynamics
clear
clc
close all
figure('position',[188,100,1300/2,550])
plot(2:9,get_aps(2:9),'Color','k','LineWidth',4)
xlabel('Time (h)')
ylabel('Concentration (μM)')
set(gca,'FontSize',28)
xlim([2 9])
% legend('WT')
f = gcf;
exportgraphics(f,'APDynamics_WT.tif','Resolution',600)