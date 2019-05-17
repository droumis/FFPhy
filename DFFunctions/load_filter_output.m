function data = load_filter_output(filtOutputDirectory, filename, animals, varargin)

filetail = '.mat';
if ~isempty(varargin)
    assign(varargin{:})
end
% if ~isempty(animal)
%     animal = ['_' animal];
% end
for an = 1:length(animals)
    animal = ['_' animals{an}];
    data{an} = load(sprintf('%s/%s%s%s',filtOutputDirectory, filename, animal, filetail));
    fprintf('loaded: %s/%s%s%s \n',filtOutputDirectory, filename, animal, ...
        filetail)
end