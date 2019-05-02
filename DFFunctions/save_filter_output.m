function save_filter_output(F, filtOutputDirectory, filename, varargin)

filetail = '';
if ~isempty(varargin)
    assign(varargin{:})
end

if ~isdir(filtOutputDirectory);
    mkdir(filtOutputDirectory);
end
try
    save(sprintf('%s/%s%s',filtOutputDirectory, filename, filetail), 'F');
catch
    fprintf('saving as v7.3 instead \n')
    save(sprintf('%s/%s%s',filtOutputDirectory, filename, filetail), 'F', '-v7.3');
end
disp(sprintf('filteroutput saved to %s/%s%s', filtOutputDirectory, filename, filetail))
end