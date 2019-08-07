function save_data(F, filtOutputDirectory, filename, varargin)
% F is a struct with at least an 'animal' field

filetail = '';
per_animal = 1;
filetailpos = 0; % if 1, at end, if 0, at beginning. at end has leading underscore
if ~isempty(varargin)
    assign(varargin{:});
end
if ~isempty(filetail)
    filename = [filename, filetail];
end
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
        save_F(F(an), filtOutputDirectory, filename, animals{an}, 'filetailpos', filetailpos);
    end
else
    save_F(F, filtOutputDirectory, filename, strjoin(animals,'-'),'filetailpos', filetailpos);
end
fprintf('saving data took  %.03f sec \n', toc);
end

function save_F(F, filtOutputDirectory, filename, filetail, varargin)
filetailpos = 0;
if ~isempty(varargin)
    assign(varargin{:})
end
if ~isdir(filtOutputDirectory)
    mkdir(filtOutputDirectory);
end
if filetailpos % filterframework results format
    savepath = sprintf('%s/%s_%s.mat',filtOutputDirectory, filename, filetail);
else % filterframework data format
    savepath = sprintf('%s/%s%s.mat',filtOutputDirectory, filetail, filename);
end
save(savepath, 'F', '-v7.3');
fprintf('saved to %s \n', savepath)
end