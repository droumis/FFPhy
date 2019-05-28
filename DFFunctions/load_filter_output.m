function data = load_filter_output(filtOutputDirectory, filename, animals, varargin)
% load filter framework output per animal
% useful to make paths with make_paths, then:
% load_filter_output(paths.filtOutputDirectory, paths.filenamesave, animals)
% specify additional filename string with varargin 'filetai'

filetail = '.mat';
if ~isempty(varargin)
    assign(varargin{:})
end

for an = 1:length(animals)
    animal = ['_' animals{an}];
    data{an} = load(sprintf('%s/%s%s%s',filtOutputDirectory, filename, animal, filetail));
    fprintf('loaded: %s/%s%s%s \n',filtOutputDirectory, filename, animal, ...
        filetail)
end