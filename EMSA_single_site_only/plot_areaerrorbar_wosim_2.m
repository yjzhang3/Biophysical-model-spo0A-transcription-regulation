function plot_areaerrorbar_wosim_2(TF_conc,real_data,real_data_std,mycolors)
% data includes all strains at all times, plot patch for each strain so
% have to use a for loop
% data dimension is 8 columns and 9 time points

TF_new = [TF_conc,flip(TF_conc)]; % TF_conc must be 1*N, real data must be N*1
for i = 1:length(real_data(1,:))
    patch = fill(TF_new,[real_data(:,i)+real_data_std(:,i);flip(real_data(:,i)-real_data_std(:,i),1)],mycolors(i,:),'EdgeColor','none','FaceAlpha','0.4');
%     hold on
%     h = errorbar(TF_conc,real_data(:,i),real_data_std(:,i),'LineWidth',2,'Color',mycolors(i,:));
%     alpha = 0.25;   
%     % Set transparency (undocumented)
%     set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
%     set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
end