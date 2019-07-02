function [out] = loadAS(animal, nt, waveSet, astype)

pathdef = animaldef(lower('Demetris'));
dirstr = sprintf('%s/analyticSignal/%s', pathdef{2}, animal);
ASsavestr = sprintf('%s/nt%02d_waveSet-%s_AS.mat',dirstr, nt, waveSet);
loadstr = sprintf('%s/nt%02d_waveSet-%s_%s.mat',dirstr, nt, waveSet, astype);
out = load(loadstr);
fprintf('loading ++++++++++ %s \n',loadstr) 

end