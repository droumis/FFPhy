function [out] = loadAS(animal, nt, waveSet, astype)

animalinfo = animaldef(lower(animal));
animalID = animalinfo{1,3}; %use anim prefix for name
FFanimdir =  sprintf('%s',animalinfo{1,2});
dtypedir = 'analyticSignal';
dirstr = sprintf('%s%s/', FFanimdir, dtypedir);

loadstr = sprintf('%s/nt%02d_waveSet-%s_%s.mat',dirstr, nt, waveSet, astype);
out = load(loadstr);
disp(sprintf('loading ++++++++++ %s',loadstr)) 

end