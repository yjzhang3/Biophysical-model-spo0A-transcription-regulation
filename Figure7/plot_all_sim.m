function plot_all_sim(TF_conc,sim_data,mycolors)

for i = 1:length(sim_data(1,:))
plot(TF_conc,sim_data(:,i),'LineWidth',7.4,'LineStyle','--','Color',mycolors(i,:))
end