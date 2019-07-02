function save_data(F, filtOutputDirectory, filename, varargin)
% F is a struct with at least an 'animal' field

filetail = '';
per_animal = 1;
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
    animals = cellfun(@(x) x{1}, {F.animal}, 'un', 0);
end

if per_animal
    for an = 1:length(animals)
        save_F(F(an), filtOutputDirectory, filename, animals{an});
    end
else
    save_F(F, filtOutputDirectory, filename, strjoin(animals,'-'))
end
end

function save_F(F, filtOutputDirectory, filename, filetail)
if ~isdir(filtOutputDirectory)
    mkdir(filtOutputDirectory);
end
save(sprintf('%s/%s_%s.mat',filtOutputDirectory, filename, filetail), 'F', '-v7.3');
fprintf('saved to %s/%s_%s \n', filtOutputDirectory, filename, filetail)
end