function Main(ID,maxtime)
tic
group_array_s.g1 = 1:8;
group_array_v.g1 = 1:8;
group_array_sv.g1 = 1:8;
strain_type = 1:8; % type of strains included to optimize the model
file = 'purethermo_Sept_2024_extrarun_highera.txt';
[~,~] = estimate_energy(...
    group_array_s,group_array_v,group_array_sv,file,ID,maxtime,strain_type);
toc
end