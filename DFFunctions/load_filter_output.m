function F = load_filter_output(filtOutputDirectory, filename, varargin)

filetail = '';
if ~isempty(varargin)
    assign(varargin{:})
end
load(sprintf('%s/%s',filtOutputDirectory, filename, filetail));
fprintf('loaded: %s/%s \n',filtOutputDirectory, filename, filetail))