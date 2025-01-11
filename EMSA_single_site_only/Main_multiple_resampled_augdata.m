%% augment to obtain 1000 datasets, with the same number of replicate experimental measurements for each strain
% then fit the model to each augmented dataset.
clear
clc
%%
n = 4;
range1 = [-5,5]; % optimization bound
fn_Khalf = sprintf('aug_EMSA_n%1d_nonorm_singlesitecurvefit_inunitKhalf.txt',n); 
fn_G = sprintf('aug_EMSA_n%1d_nonorm_singlesitecurvefit_inunitG.txt',n); 
N_sample = 1000;
%% augment to obtain 1000 datasets and 1000 sets of optimized Kd
for g = 1:N_sample
    ID = round(abs(g*9-67+g^0.2+g^3.1));
    fit_augdataset(n,range1,fn_G,fn_Khalf,ID)
end

%% read all optimized parameter sets in G and save to be used in biophysical models
%%%%!!!! can directly load this section while skipping the last section to
%%%%save time with pre-saved text files!!!%%%%
A = readmatrix(fn_G);
pars_range = A(:,1:end-1); % 1000 rows (1000 resampled datasets), 3 columns (3 parameters)
save('aug_energy.mat','pars_range')

%% read all optimized parameter sets in Khalf and compute 95% confidence interval
% find the 25th and 975th by ranking 1000 parameters
A = readmatrix(fn_Khalf);
pars_range = A(:,1:end-1); % 1000 rows (1000 resampled datasets), 3 columns (3 parameters)
m = mean(pars_range,1);
Khalf_1 = sort(pars_range(:,1),'ascend');
Khalf_2 = sort(pars_range(:,2),'ascend');
Khalf_3 = sort(pars_range(:,3),'ascend');

Khalf_1_lw = Khalf_1(25);
Khalf_1_hi = Khalf_1(975);

Khalf_2_lw = Khalf_2(25);
Khalf_2_hi = Khalf_2(975);

Khalf_3_lw = Khalf_3(25);
Khalf_3_hi = Khalf_3(975);

neg = fliplr(m)-[Khalf_3_lw,Khalf_2_lw,Khalf_1_lw]; % 0A3, 0A2, 0A1
pos = [Khalf_3_hi,Khalf_2_hi,Khalf_1_hi]-fliplr(m);

%% plot with confidence interval

figure('position',[188,100,1300/2,550])
% subplot(2,1,1)
colororder([0.5,0.5,0.5])
X = categorical({'3','2','1'});
X = reordercats(X,{'3','2','1'});
Y = fliplr(m);
b = barh(X,Y);

hold on

er = errorbar(Y, X, neg, pos, '.r','horizontal');
er.Color = [0 0 0];   
er.LineWidth = 3;

ylabel('0A box')
xlabel({'K_{half} (μM)'})
set(gca,'FontSize',28)
xlim([0 3.3])
saveas(gcf,sprintf('aug_EMSA_n%1d_singlesitereportKhalf.fig',n))
f = gcf;
exportgraphics(f,sprintf('aug_EMSA_n%1d_singlesitereportKhalf.tif',n),'Resolution',300)



