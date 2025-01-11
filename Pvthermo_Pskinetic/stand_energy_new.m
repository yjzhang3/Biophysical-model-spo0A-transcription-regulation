function G0 = stand_energy_new(config,energyi)
% input:
% energyi: where first n elements are
% binding energy of each site, and the rest are interaction energy for
% every possible pair. The pair should follow the order, 12,13,14..1n,
% 23,24,25,..2n,
% and every possible triplet: 

nbd = length(config);
energy_per_site = energyi(1:nbd);
int_energy_per_pair = energyi(nbd+1:nbd+nchoosek(nbd,2));
int_energy_per_triplet = energyi(nbd+nchoosek(nbd,2)+1:nbd+nchoosek(nbd,2)+nchoosek(nbd,3));
int_energy_per_q = energyi(nbd+nchoosek(nbd,2)+nchoosek(nbd,3)+1:nbd+nchoosek(nbd,2)+nchoosek(nbd,3)+nchoosek(nbd,4));
%% adding energy of each site
G0 = sum(config.*energy_per_site, 'all' ); 
% if listi(i) = 0, meaning nothing is bound, no energy contribution from
% there (times 0 is 0), thus we need element multiplication

%% adding energy due to TF-TF or TF-mRNA interactions
% for n molecules (including mRNA), there could be n(n-1)/2 possible
% pair of interactions. 
% if a pair doesn't interact, then the energy is 0

counter = 0; % count the number of unique pairs
int_e = 0;
for ii = 1:nbd-1

    for jj = ii+1:nbd
        counter = counter+1;
        int_e = int_e + config(ii)*config(jj)*int_energy_per_pair(counter);
    end

end
G0 = G0+int_e;

counter = 0;
int_e_triple = 0;
for k1 = 1:nbd-1
    for k2 = k1+1:nbd
        for k3 = k2+1:nbd
            counter = counter+1;
            int_e_triple = int_e_triple+config(k1)*config(k2)*config(k3)*int_energy_per_triplet(counter);
%             disp([k1,k2,k3])
        end
    end
end 

G0 = G0+int_e_triple;

counter = 0;
int_e_q = 0;
for k1 = 1:nbd-1
    for k2 = k1+1:nbd
        for k3 = k2+1:nbd
            for k4 = k3+1:nbd
            counter = counter+1;
            int_e_q= int_e_q+config(k1)*config(k2)*config(k3)*config(4)*int_energy_per_q(counter);
%             disp([k1,k2,k3])
            end
        end
    end
end 
G0 = G0+int_e_q;
end