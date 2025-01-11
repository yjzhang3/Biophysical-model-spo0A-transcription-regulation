function create_parameter_file(file, fit, score,type)
% create_parameter_file: write parameter sets to file using comma-separated values
%   create_parameter_file(file, ensemble, scores)
% inputs:
%   file = output file path (csv)
%   fit = fitted parameters (vector)
%   scores = array of particle objective function scores
% outputs:
%   none

% write parameter set to file
% file_id = fopen(file, 'w');
file_id = fopen(file,type);

for ii=1:length(fit)
    fprintf(file_id, '%4f,', fit(ii));
end

fprintf(file_id, '%4f\n', score);

% close the file
fclose(file_id);
end
