function [x,fval] = fit_data(...
    nbd,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv,...
    nvars,lb,ub,file,ID)
% input: real data for all three datasets, time dependent TF concentration, 
% number of binding
% sites, time-dependent RNAPA concentration, mutation matrix

% can also specify lower and upper bound

% Energyi and vmax and RNApH are the parameter to find for

% this is written for slurm so only iteration per function

%% parameters
rng(ID,"twister")

fun = @(p) objective_function(...
    nbd,p,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv);


%% particle swarm
% options = optimoptions(@particleswarm,'Display','iter','InitialSwarmMatrix',init_point,'OutputFcn',@outfun);
options = optimoptions(@particleswarm,'Display','iter','OutputFcn',@outfun);
[x,fval,~,~] = particleswarm(fun,nvars,lb,ub,options);
create_parameter_file(file, x, fval,'a+');


function stop = outfun(optimValues,state)
stop = false; % This function does not stop the solver
persistent psw_position_file
persistent hist
stop = false;
switch state
    case 'init'
        hist = [0,optimValues.bestx,optimValues.bestfval];
    case 'iter'
        nextline = [optimValues.iteration,optimValues.bestx,optimValues.bestfval];
%         disp(nextline)
        hist = [hist;nextline];
        psw_position_file = 'psw_position.txt';
        create_parameter_file(psw_position_file, optimValues.bestx, optimValues.bestfval,'w');
    case 'done'
        assignin('base','hist',hist);
        assignin('base','psw_position_file',psw_position_file);
end
end

end
