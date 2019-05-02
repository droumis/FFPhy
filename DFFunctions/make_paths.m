function p = make_paths(Fp)

pathdef = animaldef(lower('Demetris'));

p.filtOutputDirectory = sprintf('%s%s/', pathdef{2}, Fp.filtfunction);
p.figdirectory = sprintf('%s%s/', pathdef{4}, Fp.filtfunction);
p.resultsDirectory = sprintf('%s%s/', pathdef{3}, Fp.filtfunction);

p.filenamesave = sprintf('%s_%ss_%s', strjoin(Fp.animals,'-'), ...
    Fp.filtfunction, strjoin(Fp.epochEnvironment,'-'));
% p.filename = sprintf('%s_%s.mat', p.filenamesave, Fp.filtfunction);

end