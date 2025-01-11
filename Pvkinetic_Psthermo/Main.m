% Ps is thermo, Pv is kinetic
function Main(ID,maxtime)
group_array_s.g1 = 1:8;

% non-multiplicative vmax grouping
group_array_v.g1 = 1; %000
group_array_v.g2 = 2; %001
group_array_v.g3 = 3; 
group_array_v.g4 = 4; 
group_array_v.g5 = 5; 
group_array_v.g6 = 6; 
group_array_v.g7 = 7; 
group_array_v.g8 = 8; 

group_array_sv.g1 = 1; 
group_array_sv.g2 = 2; 
group_array_sv.g3 = 3; 
group_array_sv.g4 = 4; 
group_array_sv.g5 = 5; 
group_array_sv.g6 = 6; 
group_array_sv.g7 = 7; 
group_array_sv.g8 = 8; 

file = 'pvkinetic_psthermo_Sept_2024_extrarun_highera.txt';
[~,~] = estimate_energy(...
    group_array_s,group_array_v,group_array_sv,file,ID,maxtime);
end