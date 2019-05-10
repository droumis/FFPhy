function F = load_filter_output(filtOutputDirectory, filename, varargin)

filetail = '.mat';
if ~isempty(varargin)
    assign(varargin{:})
end
load(sprintf('%s/%s%s',filtOutputDirectory, filename, filetail));
fprintf('loaded: %s/%s%s \n',filtOutputDirectory, filename, filetail)