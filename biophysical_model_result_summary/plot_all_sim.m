function plot_all_sim(sim_data,mycolors)

for i = 1:length(sim_data(1,:))
plot(2:9,sim_data(:,i),'LineWidth',5.5,'LineStyle','--','Color',mycolors(i,:))
hold on
end