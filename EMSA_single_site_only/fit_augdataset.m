function fit_augdataset(n,range1,fn_G,fn_Khalf,ID)

%% load experimental data
load('raw_EMSA_Data_Dec_2023.mat')
% there is a total of 20*11 data points. Don't generate much more than this
% number in each augmented datasets
n_exp_per_strain = [4,3,3,3,3,2,2]; % number of experiments collected for each strain

%% allocate space for the new dataset
A1_present_aug = zeros(11,n_exp_per_strain(1)); % every n_exp_per_strain(1) is a newly generated spo0A-dependent data, there are N generations
A2_present_aug = zeros(11,n_exp_per_strain(2));
A3_present_aug = zeros(11,n_exp_per_strain(3));

%% generate one new dataset
A1_present_aug(:,1:n_exp_per_strain(1)) = normal_randi_EMSA_per_strain(A1_present,ID,n_exp_per_strain(1));
A2_present_aug(:,1:n_exp_per_strain(2)) = normal_randi_EMSA_per_strain(A2_present,ID,n_exp_per_strain(2));
A3_present_aug(:,1:n_exp_per_strain(3)) = normal_randi_EMSA_per_strain(A3_present,ID,n_exp_per_strain(3));


%% construct real and sim dataset to be fed into objective function
real_data = [A1_present_aug,...
    A2_present_aug,...
    A3_present_aug];

%% optimize
TF_conc = 0:0.2:2;

fun1 = @(G) objective_G_together_strain_weight(n,G,TF_conc,real_data,n_exp_per_strain);
nvars = 3;
lb = zeros(1,3)+range1(1);
ub = zeros(1,3)+range1(2);

x_all = zeros(20,nvars);
err_all = zeros(20,1);
for iu = 1:20
[x,err] = particleswarm(fun1,nvars,lb,ub);
x_all(iu,:) = x;
err_all(iu) = err;
end

[minerr,I] = min(err_all); % minimum error and where does it occur
create_parameter_file(fn_G, x_all(I,:), minerr,'a+'); % store in unit of free energy G

Khal_uM = exp(x_all(I,:)/4);
create_parameter_file(fn_Khalf, Khal_uM, minerr,'a+'); % store in unit of Khalf (um)


end


