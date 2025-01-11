function G0 = stand_energy_linear_new(config,energyi)
% input:
% energyi: where first n elements are
% binding energy of each site, and the rest are interaction energy for
% every possible pair. The pair should follow the order, 12,13,14..1n,
% 23,24,25,..2n,
% and every possible triplet.
% note that in this file, each parameter is exp(-G) (the negative sign is
% included!!)

nbd = length(config);
energy_per_site = energyi(1:nbd);
int_energy_per_pair = energyi(nbd+1:nbd+nchoosek(nbd,2));
int_energy_per_triplet = energyi(nbd+nchoosek(nbd,2)+1:nbd+nchoosek(nbd,2)+nchoosek(nbd,3));
% int_energy_per_q = energyi(nbd+nchoosek(nbd,2)+nchoosek(nbd,3)+1:nbd+nchoosek(nbd,2)+nchoosek(nbd,3)+nchoosek(nbd,4));

%% multiplying energy of each site
G0 = 1;
for i = 1:length(config)
    if config(i) == 1
        G0 = G0*energy_per_site(i);
    end
end
% if listi(i) = 0, meaning nothing is bound, no energy contribution from
% there (times 0 is 0), thus we need element multiplication

%% multiplying energy due to TF-TF or TF-mRNA interactions

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

%% multiplying energy due to tertiary interactions
counter = 0;
int_e_triple = 1;
for k1 = 1:nbd-1
    for k2 = k1+1:nbd
        for k3 = k2+1:nbd
            counter = counter+1;
            if (config(k1) == 1) && (config(k2) == 1) && (config(k3) == 1)
                int_e_triple = int_e_triple*int_energy_per_triplet(counter);
            end
        end
    end
end 

G0 = G0*int_e_triple;
% 
% %% multiply energy due to quarternery interactons
% counter = 0;
% int_e_q = 1;
% for k1 = 1:nbd-1
%     for k2 = k1+1:nbd
%         for k3 = k2+1:nbd
%             for k4 = k3+1:nbd
%                 counter = counter+1;
%                 if (config(k1) == 1) && (config(k2) == 1) && (config(k3) == 1) && (config(k4) == 1)
%                     int_e_q = int_e_q*int_energy_per_q(counter);
%                 end
%             end
%         end
%     end
% end 
% 
% G0 = G0*int_e_q;

end