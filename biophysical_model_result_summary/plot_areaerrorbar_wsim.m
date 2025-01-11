function plot_areaerrorbar_wsim(real_data,real_data_std,sim_data,mycolors)
% data includes all strains at all times, plot patch for each strain so
% have to use a for loop
% data dimension is 8 columns and 9 time points

for i = 1:length(real_data(1,:))
    patch = fill([2 3 4 5 6 7 8 9 9 8 7 6 5 4 3 2], [real_data(:,i)+real_data_std(:,i);flip(real_data(:,i)-real_data_std(:,i),1)],mycolors(i,:),'EdgeColor','none','FaceAlpha','0.6');
    hold on
    h = errorbar(2:9,real_data(:,i),real_data_std(:,i),'.','LineStyle','none','Color',mycolors(i,:));
    % Set transparency level (0:1)
    alpha = 0.5;   
    % Set transparency (undocumented)
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
    set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])

    hold on
    c = [mycolors(i,:),1];
    plot(2:9,sim_data(:,i),'LineWidth',2,'LineStyle','--','Color',c)
end