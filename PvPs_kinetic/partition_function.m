function Z = partition_function(nbd,energyi,mut,TF_conc)
% calculate partition function of the configurations as if promoter is not
% there.

nbd_0A = nbd; % I've already excluded the promoter in the higher file (transciption_rate_new)!

bins = dec2bin(0:(2^nbd_0A-1), nbd_0A) - '0'; % all possible configurations when promoter is not there
% 2^nbd rows representing n possible configurations, nbd columns represent # of 0A
% boxes

% consider mutated sites
for ss = 1:length(mut)
    if mut(ss) == 0
        on_ind = find(bins(:,ss)~=1);
        bins = bins(on_ind,:);
    end
end

Z = 0;
for bb = 1:length(bins(:,1))
    config = bins(bb,:);
    BF = Boltzmann_factor_per_config(config, energyi,TF_conc);
%     disp(BF)
    Z = Z + BF;
end

% [1,1,1] all working
% [0,0,0] all mutated
end
