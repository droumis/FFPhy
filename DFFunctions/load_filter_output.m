function data = load_filter_output(filtOutputDirectory, filename, animals, varargin)
% d = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filename, Fp.animals);
% load filter framework results

% specify additional filename string with varargin 'filetail'

% Author: Demetris Roumis June 2019

filetail = '';
if ~isempty(varargin)
    assign(varargin{:})
end

for an = 1:length(animals)
    tmp = load(sprintf('%s/%s%s_%s.mat',filtOutputDirectory, filename, filetail, ...
        animals{an}));    
    fprintf('loaded: %s/%s%s_%s.mat \n',filtOutputDirectory, filename, filetail, ...
        animals{an})
    data{an} = tmp.F;
end
data = [data{:}]'; % merge the cell array of structs to match filter output form