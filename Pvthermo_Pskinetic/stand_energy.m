function G0 = stand_energy(config,energyi)
% input:
% energyi: where first n elements are
% binding energy of each site, and the rest are interaction energy for
% every possible pair. The pair should follow the order, 12,13,14..1n,
% 23,24,25,..2n,

nbd = length(config);
energy_per_site = energyi(1:nbd);
int_energy_per_pair = energyi(nbd+1:end);


%% adding energy of each site
G0 = sum(config.*energy_per_site, 'all' ); 
% if listi(i) = 0, meaning nothing is bound, no energy contribution from
% there (times 0 is 0), thus we need element multiplication

%% adding energy due to TF-TF or TF-mRNA interactions
% for n molecules (including mRNA), there could be n(n-1)/2 possible
% pair of interactions. 
% if a pair doesn't interact, then the energy is 0

count = 0; % count the number of unique pairs
int_e = 0;
for ii = 1:nbd-1

    for jj = ii+1:nbd
        count = count+1;
        int_e = int_e + config(ii)*config(jj)*int_energy_per_pair(count);
    end

end
G0 = G0+int_e;

end