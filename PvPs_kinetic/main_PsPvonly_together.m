function main_PsPvonly_together(ID,maxtime,filename)
%%
group_array_s.g1 = 1; 
group_array_s.g2 = 2; 
group_array_s.g3 = 3; 
group_array_s.g4 = 4; 
group_array_s.g5 = 5; 
group_array_s.g6 = 6; 
group_array_s.g7 = 7; 
group_array_s.g8 = 8; 

%% Pv grouping
group_array_v.g1 = 1; 
group_array_v.g2 = 2; 
group_array_v.g3 = 3; 
group_array_v.g4 = 4; 
group_array_v.g5 = 5; 
group_array_v.g6 = 6; 
group_array_v.g7 = 7; 
group_array_v.g8 = 8; 

%% PvPs grouping
group_array_sv.g1 = 1; 
group_array_sv.g2 = 2; 
group_array_sv.g3 = 3; 
group_array_sv.g4 = 4; 
group_array_sv.g5 = 5; 
group_array_sv.g6 = 6; 
group_array_sv.g7 = 7; 
group_array_sv.g8 = 8; 


%% finish set up

[~,~] = estimate_energy_together(group_array_s,group_array_v,group_array_sv,filename,ID,maxtime);

end
