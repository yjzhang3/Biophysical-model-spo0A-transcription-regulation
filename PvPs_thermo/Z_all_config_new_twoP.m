function Zall = Z_all_config_new_twoP(nbd,energyi,mut,TF_conc,RNApH_conc,RNApA_conc)
% denoinator should include all possible configurations (even if RNAp is
% not bound) but of course mutation would make some configurations have a
% boltzmann factor of 0, so it's equivalent to not including them

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
    Z = Z_per_config_new_twoP(config,energyi,TF_conc,RNApH_conc,RNApA_conc);
    Zall = Zall + Z;
end


end