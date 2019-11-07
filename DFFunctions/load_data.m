function out = load_data(filtOutputDirectory, filename, animals, varargin)
% d = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filename, Fp.animals);
% load filter framework results

% specify additional filename string with varargin 'filetail'
if ~isa(animals, 'cell')
    animals = {animals};
end
% Author: Demetris Roumis June 2019
filtfunction = ' ';
filetail = '';
animpos = 1; % animal name position. 0 at beginning without underscore. 1 for end with it
if ~isempty(varargin)
    assign(varargin{:})
end
tic
fprintf('loading %s\n', animals{:})
for an = 1:length(animals)
    animal = animals{an};
    andef = animaldef(animal);
    if strcmp(filtOutputDirectory, 'filterframework')
        dirout = andef{2};
    else
        dirout = filtOutputDirectory;
    end
    switch filtfunction
        case 'behavestate'
            tmp = load(sprintf('%s%sBehaveState.mat',andef{2}, andef{3}));
            fprintf('%s%sBehaveState.mat \n',andef{2}, andef{3});
            F = tmp.BehaveState;
        otherwise
            if animpos
                tmp = load(sprintf('%s/%s%s_%s.mat',dirout , filename, filetail, animal));
                fprintf('loaded: %s/%s%s_%s.mat \n',dirout , filename, filetail, animal)
            else
                tmp = load(sprintf('%s/%s%s%s.mat', dirout, animal, filename, filetail));
                fprintf('loaded: %s/%s%s%s.mat \n',animal, dirout , filename, filetail)
            end
            F = tmp.F;
    end
    out(an) = F(1);
    out(an).animal = andef;
end
% data = [data{:}]'; % merge the cell array of structs to match filter output form
fprintf('took %.02f sec\n',toc);
end