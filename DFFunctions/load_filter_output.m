function F = load_filter_output(filtOutputDirectory, filename, varargin)

filetail = '.mat';
animal = '';
if ~isempty(varargin)
    assign(varargin{:})
end
if ~isempty(animal)
    animal = ['_' animal];
end
load(sprintf('%s/%s%s%s',filtOutputDirectory, filename, animal, filetail));
fprintf('loaded: %s/%s%s%s \n',filtOutputDirectory, filename, animal, filetail)