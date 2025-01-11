function aug_data = normal_randi_EMSA_per_strain(data,ID,n) % for a particular strain
% first column is concentration, second column is mean, third column is std
% n = generate random number n times accroding to how many experiments are
% actualy conducted for this strain

% output is 11 rows by n columns, each n is a random draw
rng(ID,'twister')
aug_data = zeros(length(data(:,1)),n);
for j = 1:n
    aug_data_j = zeros(length(data(:,1)),1);
    for i = 1:length(data(:,1))
    %     disp(data(i,2))
    %     disp(data(i,3))
    %     disp(normrnd(data(i,2),data(i,3)))
        aug_data_j(i) = normrnd(data(i,2),data(i,3));
    end
    aug_data(:,j) = aug_data_j;
end
    
