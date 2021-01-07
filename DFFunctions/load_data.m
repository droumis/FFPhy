function out = load_data(source, filename, animals, varargin)
% load filter framework results
% source can be full path OR one of a special set: 'filterframework'
% >> d = (Fp.paths.filtOutputDirectory, Fp.paths.filename, Fp.animals);
% >> dio = (Fp.paths.filtOutputDirectory, Fp.paths.filename, Fp.animals);
% 

% specify additional filename string with varargin 'filetail'
if ~isa(animals, 'cell')
    animals = {animals};
end

pconf = paramconfig;
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
    if strcmp(source, 'filterframework')
        dirout = andef{2};
%     elseif strcmp(source, 'results') % now in Fp.paths.resultsDirectory
%         dirout = [pconf.andef{3} '/' filename, '/'];
    else
        dirout = source;
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
            try
                F = tmp.F; % filter output struct
            catch
                F = tmp; % anything else
            end
    end
    try
        out(an) = F(1);
        out(an).animal = andef;
    catch
        out = tmp; % assumes loading one animal
    end
end
% data = [data{:}]'; % merge the cell array of structs to match filter output form
fprintf('took %.02f sec\n',toc);
end