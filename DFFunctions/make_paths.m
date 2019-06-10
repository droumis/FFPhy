function p = make_paths(filtfunction, varargin)

% epochEnvironment = 'alleps';
if ~isempty(varargin)
    assign(varargin{:})
end
pathdef = animaldef(lower('Demetris'));

p.filtOutputDirectory = sprintf('%s%s/', pathdef{2}, filtfunction);
p.figdirectory = sprintf('%s%s/', pathdef{4}, filtfunction);
p.resultsDirectory = sprintf('%s%s/', pathdef{3}, filtfunction);
p.filenamesave = filtfunction;
% p.filename = sprintf('%s_%s.mat', p.filenamesave, Fp.filtfunction);

end