function save_data(F, filtOutputDirectory, filename, varargin)
% F is a struct with at least an 'animal' field
pconf = paramconfig;
filetail = '';
per_animal = 1;
animpos = 1; % if 1, at end, if 0, at beginning. at end has leading underscore
varname = '';
if ~isempty(varargin)
    assign(varargin{:});
end
% if ~isempty(filetail)
%     filename = [filename, filetail];
% end
% the animal field should either just be their name or the loaded content
% of andef (currently 6 cell array)
try
        animals = cellfun(@(x) x{3}, {F.animal}, 'un', 0);
catch
    try
        animals = cellfun(@(x) x{1}, {F.animal}, 'un', 0);
    catch
        animals = cellfun(@(x) x, {F.animal}, 'un', 0);
    end
end
tic
if per_animal
    for an = 1:length(animals)
        if strcmp(filtOutputDirectory, 'filterframework')
            andef = animaldef(animals{an});
            dirout = andef{2};
            save_F(F(an), dirout, [filename filetail], animals{an}, ...
                'animpos', animpos, 'varname', varname);
        elseif strcmp(filtOutputDirectory, 'results')
            dirout = [pconf.andef{3} '/' filename, '/'];
            save_F(F(an), dirout, [filename filetail], animals{an}, ...
                'animpos', animpos, 'varname', varname);
        else
            dirout = filtOutputDirectory;        
            save_F(F(an), dirout, [filename filetail], animals{an}, ...
                'animpos', animpos, 'varname', varname);
        end
    end
else
    save_F(F, filtOutputDirectory, [filename filetail], strjoin(animals,'-'), ...
        'animpos', animpos, 'varname', varname);
end
fprintf('saving data took  %.03f sec \n', toc);
end

function save_F(F, filtOutputDirectory, filename, animal, varargin)
animpos = 1;
varname = '';
if ~isempty(varargin)
    assign(varargin{:});
end
if ~isdir(filtOutputDirectory)
    mkdir(filtOutputDirectory);
end
if animpos % results format
    savepath = sprintf('%s/%s_%s.mat',filtOutputDirectory, filename, animal);
else % data format
    savepath = sprintf('%s/%s%s.mat',filtOutputDirectory, animal, filename);
end
if any(varname)
    eval([varname, ' = F']);
    save(savepath, varname, '-v7.3');
else
    save(savepath, 'F', '-v7.3');
end
fprintf('saved to %s \n', savepath)
end