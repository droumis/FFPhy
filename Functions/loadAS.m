function [out] = loadAS(animal, nt, waveSet, lfptype)

pathdef = animaldef(lower('Demetris'));
dirstr = sprintf('%s/analyticSignal/%s', pathdef{2}, animal);
% loadstr = sprintf('%s/nt%02d_waveSet-%s_%s.mat',dirstr, nt, waveSet, astype);
loadstr = sprintf('%s/AnalyticSignal_nt%02d_waveSet-%s_%s.mat', dirstr, ...
    nt, waveSet, lfptype);

out = load(loadstr);
fprintf('loading ++++++++++ %s \n',loadstr) 

end