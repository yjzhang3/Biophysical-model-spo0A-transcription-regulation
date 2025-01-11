% compute BIC for purely kinetic model with different number of energy
% parameters being zero

% BIC = log(𝑆𝑆𝐸)+𝑁_𝑝𝑎𝑟𝑠/𝑁_𝑜𝑏𝑠 ×log(𝑁_𝑜𝑏𝑠)
clear
hexColor=['#e60049';'#50e991';'#e6d800';'#9b19f5';'#ffa300';'#dc0ab4';'#00bfa0';'#0bb4ff'];
rgbColor = hex2rgb(hexColor);

% all parameters are 0 
fn = 'local_run_61000_121323123zero.txt';
A = readmatrix(fn);
[M_1,~] = min(A(:,end));

% 3 parameters are 0

fn = 'local_run_61000_1323123zero.txt';
A = readmatrix(fn);
[M_2,~] = min(A(:,end));

fn = 'local_run_61000_1223123zero.txt';
A = readmatrix(fn);
[M_3,~] = min(A(:,end));

fn = 'local_run_61000_1213123zero.txt';
A =readmatrix(fn);
[M_5,~] = min(A(:,end));

fn = 'local_run_61000_121323zero.txt';
A = readmatrix(fn);
[M_6,~] = min(A(:,end));

% 2 parameters are 0
fn = 'local_run_61000_23123zero.txt';
A = readmatrix(fn);
[M_7,~] = min(A(:,end));

fn = 'local_run_61000_13123zero.txt';
A = readmatrix(fn);
[M_8,~] = min(A(:,end));

fn = 'local_run_61000_1323zero.txt';
A = readmatrix(fn);
[M_9,~] = min(A(:,end));

fn = 'local_run_61000_12123zero.txt';
A = readmatrix(fn);
[M_10,~] = min(A(:,end));

fn = 'local_run_61000_1223zero.txt';
A = readmatrix(fn);
[M_11,~] = min(A(:,end));

fn = 'local_run_61000_1213zero.txt';
A = readmatrix(fn);
[M_12,~] = min(A(:,end));

% 1 parameter is 0 
fn = 'local_run_61000_123zero.txt'; 
A = readmatrix(fn);
[M_13,~] = min(A(:,end));

fn = 'local_run_61000_23zero.txt'; 
A = readmatrix(fn);
[M_14,~] = min(A(:,end));

fn = 'local_run_61000_13zero.txt';
A = readmatrix(fn);
[M_15,~] = min(A(:,end));

fn = 'local_run_61000_12zero.txt';
A = readmatrix(fn);
[M_16,~] = min(A(:,end));

% 0 parameters are zero
fn = 'purekinetic_Sept_2024_extrarun_highera.txt';
A = readmatrix(fn);
[M_Z,~] = min(A(:,end));

N_pars0 = 0+1+32; % 32 = 2 groups of vmax and 2 groups of R, 1 = binding energy parameters, n = free interaction parameters that exist
N_pars1 = 1+1+32;
N_pars2 = 2+1+32;
N_pars3 = 3+1+32;
N_pars4 = 4+1+32;
N_obs = 8*8*3; 

BIC_n0 = log(M_1)+N_pars0/N_obs*log(N_obs);

BIC_n1 = zeros(1,4);
BIC_n1(1) = log(M_2)+N_pars1/N_obs*log(N_obs);
BIC_n1(2) = log(M_3)+N_pars1/N_obs*log(N_obs);
BIC_n1(3) = log(M_5)+N_pars1/N_obs*log(N_obs);
BIC_n1(4) = log(M_6)+N_pars1/N_obs*log(N_obs);

BIC_n2 = zeros(1,6);
BIC_n2(1) = log(M_7)+N_pars2/N_obs*log(N_obs);
BIC_n2(2) = log(M_8)+N_pars2/N_obs*log(N_obs);
BIC_n2(3) = log(M_9)+N_pars2/N_obs*log(N_obs);
BIC_n2(4) = log(M_10)+N_pars2/N_obs*log(N_obs);
BIC_n2(5) = log(M_11)+N_pars2/N_obs*log(N_obs);
BIC_n2(6) = log(M_12)+N_pars2/N_obs*log(N_obs);

BIC_n3 = zeros(1,4);
BIC_n3(1) = log(M_13)+N_pars3/N_obs*log(N_obs);
BIC_n3(2) = log(M_14)+N_pars3/N_obs*log(N_obs);
BIC_n3(3) = log(M_15)+N_pars3/N_obs*log(N_obs);
BIC_n3(4) = log(M_16)+N_pars3/N_obs*log(N_obs);

BIC_n4 = log(M_Z)+N_pars4/N_obs*log(N_obs);

%%
% close all
% figure1 = figure("position",[229,193,958,477]); % for additivity
% colororder([rgbColor(8,:)])
% 
% 
% % model_type = {'Model unconstrained with EMSA prediction';'Model constrained with EMSA prediction'};
% % Y = [BIC_1;BIC_2];
% % mltable = table(Y,'RowNames',model_type);
% 
% 
% X = categorical({'Model partially constrained with EMSA prediction','Model completely constrained with EMSA prediction'});
% X = reordercats(X,{'Model partially constrained with EMSA prediction','Model completely constrained with EMSA prediction'});
% Y = [BIC_1,BIC_2];
% b = bar(X,Y);
% 
% annotation(figure1,'textbox',...
%     [0.299371607515657 0.6180062893081762-0.05 0.197329853862214 0.0398322851153053],...
%     'String',string(BIC_1),...
%     'LineStyle','none',...
%     'FontSize',18,...
%     'FitBoxToText','off','VerticalAlignment','middle',...
%     'HorizontalAlignment','center');
% 
% annotation(figure1,'textbox',...
%     [0.563981210855949 0.751656184486373-0.05 0.206724425887266 0.0419287211740041],...
%     'String',string(BIC_2),...
%     'LineStyle','none',...
%     'FontSize',18,...
%     'FitBoxToText','off','VerticalAlignment','middle',...
%     'HorizontalAlignment','center');
% 
% % title('Parameter and Erorr Tradeoff')
% title('BIC')
% set(gca,'FontSize',20)
% % ylim([-2 0])
% 
% saveas(gcf,'BIC.fig')
% % saveas(gcf,'BIC.tif')