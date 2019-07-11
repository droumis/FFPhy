function data = load_data(filtOutputDirectory, filename, animals, varargin)
% d = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filename, Fp.animals);
% load filter framework results

% specify additional filename string with varargin 'filetail'
if ~isa(animals, 'cell')
    animals = {animals};
end
% Author: Demetris Roumis June 2019
filtfunction = ' ';
filetail = '';
if ~isempty(varargin)
    assign(varargin{:})
end
tic
fprintf('loading %s\n', animals{:})
for an = 1:length(animals)
    andef = animaldef(animals{an});
    switch filtfunction
        case 'behavestate'
            tmp = load(sprintf('%s%sBehaveState.mat',andef{2}, andef{3}));
            fprintf('%s%sBehaveState.mat \n',andef{2}, andef{3});
            data{an} = tmp.BehaveState;
        otherwise
            tmp = load(sprintf('%s/%s%s_%s.mat',filtOutputDirectory, filename, filetail, ...
                animals{an}));
            fprintf('loaded: %s/%s%s_%s.mat \n',filtOutputDirectory, filename, filetail, ...
                animals{an})
            data{an} = tmp.F;
    end
end
data = [data{:}]'; % merge the cell array of structs to match filter output form
fprintf('took %.02f sec\n',toc);
end