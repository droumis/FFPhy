function save_filter_output(F, filtOutputDirectory, filename, varargin)

per_animal = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
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
if ~isdir(filtOutputDirectory);
    mkdir(filtOutputDirectory);
end
save(sprintf('%s/%s_%s.mat',filtOutputDirectory, filename, filetail), 'F', ...
    '-v7.3');
fprintf('filteroutput saved to %s/%s_%s \n', filtOutputDirectory, filename, ...
    filetail)
end