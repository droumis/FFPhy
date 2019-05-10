function p = make_paths(filtfunction, epochEnvironment)

pathdef = animaldef(lower('Demetris'));

p.filtOutputDirectory = sprintf('%s%s/', pathdef{2}, filtfunction);
p.figdirectory = sprintf('%s%s/', pathdef{4}, filtfunction);
p.resultsDirectory = sprintf('%s%s/', pathdef{3}, filtfunction);
try
    p.filenamesave = sprintf('%s_%s', ...
        filtfunction, strjoin(epochEnvironment,'-'));
catch
    p.filenamesave = sprintf('%s_%s', ...
        filtfunction, epochEnvironment);
end
% p.filename = sprintf('%s_%s.mat', p.filenamesave, Fp.filtfunction);

end