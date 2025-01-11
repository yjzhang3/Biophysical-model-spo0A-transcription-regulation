function G0 = stand_energy_linear(config,energyi)
% input:
% energyi: where first n elements are
% binding energy of each site, and the rest are interaction energy for
% every possible pair. The pair should follow the order, 12,13,14..1n,
% 23,24,25,..2n,
% note that in this file, each parameter is exp(-G) (the negative sign is
% included!!)
% promoter is not included!

nbd = length(config);
energy_per_site = energyi(1:nbd);
int_energy_per_pair = energyi(nbd+1:end);

%% adding energy of each site
G0 = 1;
for i = 1:length(config)
    if config(i) == 1
        G0 = G0*energy_per_site(i);
    end
end
% if listi(i) = 0, meaning nothing is bound, no energy contribution from
% there (times 0 is 0), thus we need element multiplication

%% adding energy due to TF-TF or TF-mRNA interactions
% for n molecules (including mRNA), there could be n(n-1)/2 possible
% pair of interactions. 
% if a pair doesn't interact, then the energy is 0

count = 0; % count the number of unique pairs
int_e = 1;
for ii = 1:nbd-1
    for jj = ii+1:nbd
        count = count+1;
        if (config(ii) == 1) && (config(jj) == 1)
            int_e = int_e*int_energy_per_pair(count);
        end
    end
end
G0 = G0*int_e;

end


