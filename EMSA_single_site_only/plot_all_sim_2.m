function plot_all_sim_2(x,sim_data,sim_data_std,mycolors)
% plot dashed prediction with color patch showing standard deviation
% sim_data and std must be N*1 vector

x_new = [x,flip(x)]; % x must be 1*N vector
for i = 1:length(sim_data(1,:))
patch = fill(x_new, [sim_data(:,i)+sim_data_std(:,i);flip(sim_data(:,i)-sim_data_std(:,i),1)],mycolors(i,:),'EdgeColor','none','FaceAlpha','0.3');
hold on
plot(x,sim_data(:,i),'LineWidth',5.5,'LineStyle','--','Color',mycolors(i,:))
hold on
end