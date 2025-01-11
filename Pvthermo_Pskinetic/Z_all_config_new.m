function Zall = Z_all_config_new(nbd,energyi,mut,TF_conc,RNAp_conc)

bins = dec2bin(0:(2^nbd-1), nbd) - '0'; % all possible configurations

% consider mutated sites
for ss = 1:length(mut)
    if mut(ss) == 0
        on_ind = find(bins(:,ss)~=1);
        bins = bins(on_ind,:);
    end
end

Zall = 0;
for bb = 1:length(bins(:,1))
    config = bins(bb,:);
    Z = Z_per_config_new(config, energyi,TF_conc,RNAp_conc);
    Zall = Zall + Z;
end

% Zall = Zall/2^(nbd-1-sum(mut));
% divide by 2^(number of mutations), where number of mutations is 
% nbd-1-sum(mut), becauase nbd-1 is the number of binding sites besides the
% promoter, and sum(mut) is the number of non-mutated sites

% [1,1,1] all working
% [0,0,0] all mutated
end
