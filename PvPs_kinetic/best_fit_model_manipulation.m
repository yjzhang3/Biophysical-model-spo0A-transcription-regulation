%% load all data, in this file, we plot different possibilities (constant sigma factor or 0A~P ) and compare to original prediction
clear
clc
best_fit_model_prediction
best_fit_model_prediction_constant_0AP
best_fit_model_prediction_constant_pol
%% 
clear 
clc
close all
load('sim_data_const_holoenzyme.mat') % pre-run and saved
load('sim_data_const_0AP.mat')
load('sim_data_original.mat')
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

TF_conc_t = get_aps(2:9);

%% Sim data normalized (constant sigma) and nonnormalized (original sim data) for just WT, along with normzlied 123* 
load('best_fit_model_prediction_colors.mat')
mycolors = hex2rgb(hexColor);

linewid = 7.4;
%%%%%% pv
figure('position',[188,100,1300/2*1.5,550*1.5]); % for strain fits, single row subplot(2,1,1)
for i = 1:1
    plot(2:9,sim_data_v(:,8),'LineWidth',linewid,'LineStyle','--','Color',mycolors(8,:))
    hold on
    plot(2:9,sim_data_v(:,i),'LineWidth',linewid,'LineStyle','--','Color',mycolors(i,:))
    hold on
    plot(2:9,sim_data_v_constsig(:,i),'LineWidth',linewid,'LineStyle',':','Color',mycolors(i,:))
    hold on 
%     plot(2:9,sim_data_v_constTF(:,i),'LineWidth',linewid,'LineStyle','--','Marker',"x",'MarkerSize',27,'MarkerFaceColor',mycolors(i,:),'Color',mycolors(i,:))
%     hold on
    
end
ylim([0 1460])
xlim([2 9])
ax = gca;
box(ax,'on')
% box(f1,'on');
xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
[~, hobj, ~, ~] = legend({['none'],['123' newline ''],['123, constant {RNAP-\sigma^{A}}' newline '']},'Location','northwest','NumColumns',3,'FontSize',32,'Location','Northwest');
% 'Position',[0.281709401709402,0.664040411229085,0.487692297697068,0.247878780690107]
% ,['123, constant 0A~P' newline '']
% [~, hobj, ~, ~] = legend({['Dynamic {\sigma^{A}} (no 0A boxes)'],['Dynamic {\sigma^{A}}, constant 0A~P' newline ''],['Dynamic {\sigma^{A}}, dynamic 0A~P' newline ''],['Constant {\sigma^{A}}, dynamic 0A~P' newline '']},'Location','northwest','FontSize',28)
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.9);
set(gca,'FontSize',38)

f = gcf;
exportgraphics(f,'Constant_sig_WT_Pv.tif','Resolution',600)
saveas(f,'FigB_Constant_sig_WT_Pv.svg')


%%
%%%%%%%%% ps
figure('position',[188,100,1300/2*1.5,550*1.5]); % for strain fits, single row subplot(2,1,1)
% subplot3=subplot(1,3,3)
for i = 1:1
    plot(2:9,sim_data_s(:,8),'LineWidth',linewid,'LineStyle','--','Color',mycolors(8,:))    
    hold on
        plot(2:9,sim_data_s(:,i),'LineWidth',linewid,'LineStyle','--','Color',mycolors(i,:))
    hold on
    plot(2:9,sim_data_s_constsig(:,i),'LineWidth',linewid,'LineStyle',':','Color',mycolors(i,:))
    hold on 

%     plot(2:9,sim_data_s_constTF(:,i),'LineWidth',linewid,'LineStyle','--','Marker',"x",'MarkerSize',27,'MarkerFaceColor',mycolors(i,:),'Color',mycolors(i,:))
%     hold on

end
ax = gca;
box(ax,'on')
ylim([0 1460])
xlim([2 9])
xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
[~, hobj, ~, ~] = legend({['none'],['123' newline ''],['123, constant RNAP-{\sigma^{H}}' newline '']},'Location','northwest','NumColumns',3,'FontSize',32,'Location','northwest');
% 'Position',[0.281709401709402,0.664040411229085,0.487692297697068,0.247878780690107]
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.9);
set(gca,'FontSize',38)
f = gcf;
exportgraphics(f,'Constant_sig_WT_Ps.tif','Resolution',600)
saveas(f,'FigC_Constant_sig_WT_Ps.svg')

%%
%%%%%%%%%% pspv
figure('position',[188,100,1300/2*1.5,550*1.5]); % for strain fits, single row subplot(2,1,1)
% subplot1=subplot(1,3,1)
for i = 1:1
    plot(2:9,sim_data_sv(:,8),'LineWidth',linewid,'LineStyle','--','Color',mycolors(8,:))
    hold on
        plot(2:9,sim_data_sv(:,i),'LineWidth',linewid,'LineStyle','--','Color',mycolors(i,:))
    hold on
    plot(2:9,sim_data_sv_constsig(:,i),'LineWidth',linewid,'LineStyle',':','Color',mycolors(i,:))
    hold on



end
ylim([0 1460])
xlim([2 9])
[~, hobj, ~, ~] = legend({['none' newline ''],['123' newline ''],['123, constant RNAP-{\sigma^{A/H}}' newline '']},'Location','northwest','NumColumns',3,'FontSize',32,'Location','northwest');
% 'Position',[0.283376069552878,0.631767685003954,0.546666654317807,0.278181808897943]
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2.9);

xlabel('Time (h)')
ylabel('v^{effective} (a.u.)')
% title('PvPs strain')
% box(f3,'on');
set(gca,'FontSize',38)
f = gcf;
exportgraphics(f,'Constant_sig_WT_PvPs.tif','Resolution',600)
saveas(f,'FigA_Constant_sig_WT_PvPs.svg')