clear
clc

hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);


%% read data
kfn = 'purekinetic_Sept_2024_extrarun_highera.txt';
A = readmatrix(kfn);
km = mean(A(:,end));
kstd = std(A(:,end));

tfn = 'purethermo_Sept_2024_extrarun_highera.txt';
A = readmatrix(tfn);
tm = mean(A(:,end));
tstd = std(A(:,end));

ktfn = 'pvkinetic_psthermo_Sept_2024_extrarun_highera.txt';
A = readmatrix(ktfn);
ktm = mean(A(:,end));
ktstd = std(A(:,end));

tkfn = 'pvthermo_pskinetic_Sept_2024_extrarun_highera.txt';
A = readmatrix(tkfn);
tkm = mean(A(:,end));
tkstd = std(A(:,end));

%%
figure("position",[680,678,560*2,420]); % for additivity
colororder([0.5,0.5,0.5])
X = categorical({'Purely Kinetic','Purely Thermodynamic','Thermodynamic{\it Pv}, Kinetic{\it Ps}','Kinetic{\it Pv}, Thermodynamic{\it Ps}'});
X = reordercats(X,{'Purely Kinetic','Purely Thermodynamic','Thermodynamic{\it Pv}, Kinetic{\it Ps}','Kinetic{\it Pv}, Thermodynamic{\it Ps}'});
Y = [km,tm,ktm,tkm];
b = barh(X,Y);

hold on

er = errorbar(Y, X, [kstd,tstd,ktstd,tkstd], '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

% title('Model Fitting Error')
xlabel('Optimization Error')
set(gca,'FontSize',22)
xlim([0 0.5])

saveas(gcf,'model_summary.fig')
f = gcf;
exportgraphics(f,'model_summary.tif','Resolution',600)

