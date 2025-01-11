function plot_areaerrorbar_wosim(TF_conc,real_data,real_data_std,mycolors)
% data includes all strains at all times, plot patch for each strain so
% have to use a for loop
% data dimension is 8 columns and 9 time points

TF_new = [TF_conc;flip(TF_conc)];
for i = 1:length(real_data(1,:))
%     patch = fill([2 3 4 5 6 7 8 9 9 8 7 6 5 4 3 2], [real_data(:,i)+real_data_std(:,i);flip(real_data(:,i)-real_data_std(:,i),1)],mycolors(i,:),'EdgeColor','none','FaceAlpha','0.1');
%     hold on
%     h =
%     errorbar(TF_conc,real_data(:,i),real_data_std(:,i),'Marker','o','MarkerSize',7,'MarkerFaceColor',mycolors(i,:),'LineStyle','none','LineWidth',2.78,'Color',mycolors(i,:));
%     % defulat above
    h = errorbar(TF_conc,real_data(:,i),real_data_std(:,i),'Marker','o','MarkerSize',13,'MarkerFaceColor',mycolors(i,:),'LineStyle','none','LineWidth',5,'Color',mycolors(i,:));
    hold on
%     alpha = 0.5;   
%     % Set transparency (undocumented)
%     set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
%     set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
%     plot(2:9,real_data(:,i),'LineStyle','none','Color',mycolors(i,:))
end