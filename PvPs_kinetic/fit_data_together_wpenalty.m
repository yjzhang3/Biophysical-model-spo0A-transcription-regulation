function [x,fval] = fit_data_together_wpenalty(...
    nbd,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv,real_data_Ps,real_data_Pv,real_data_PsPv, ...
    lb,ub,...
    nvars,group_array_s,group_array_v,group_array_sv,filename,ID,maxtime)
%% parameters
rng(ID,'twister')

fun = @(p) objective_function(...
    nbd,p,TF_conc_t,...
    mut_mat_s,mut_mat_v,mut_mat_sv, ...
    real_data_Ps,real_data_Pv,real_data_PsPv,...
    group_array_s,group_array_v,group_array_sv);

%% set up options
% options = optimoptions('particleswarm', 'SwarmSize',50,'HybridFcn',@fmincon);

%% particle swarm
options = optimoptions(@particleswarm,'Display','iter','MaxTime',maxtime,'OutputFcn',@outfun);
[x,fval,~,~] = particleswarm(fun,nvars,lb,ub,options);
create_parameter_file(filename, x, fval,'a+');
% 

%%
% fn = strcat('checkfile',string(ID),'.mat');
% if isfile(fn)
%     disp('using checkfile')
%     [x,fval,exitflag,output] = surrogateopt(fn);
% %     disp(xval)
% %     disp(exitflag)
%     
% else
%     disp('starting fresh')
%     options = optimoptions("surrogateopt","CheckpointFile",fn);
%     [x,fval,exitflag,output] = surrogateopt(fun,lb,ub,options);
% end
% create_parameter_file(filename, x, fval,'w');

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
        psw_position_file = 'psw_position_1.txt';
        create_parameter_file(psw_position_file, optimValues.bestx, optimValues.bestfval,'w');
    case 'done'
        assignin('base','hist',hist);
        assignin('base','psw_position_file',psw_position_file);
end
end

end