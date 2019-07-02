function saveAS(as,nt,ntN,wp, animal)

pathdef = animaldef(lower('Demetris'));
dirstr = sprintf('%s/analyticSignal/%s', pathdef{2}, animal);
if ~isdir(dirstr)
    mkdir(dirstr);
end

analyticsignal.analyticsignal = as;
analyticsignal.waveparams = wp;
analyticsignal.ntrodeN = ntN;
analyticsignal.ntrode = nt;

ASsavestr = sprintf('%s/nt%02d_waveSet-%s_AS.mat', dirstr, ...
    analyticsignal.ntrode, analyticsignal.waveparams.waveSet);
save(ASsavestr, 'analyticsignal', '-v7.3');
fprintf('SAVED ANALYTIC SIGNAL RESULTS ++++++++++ %s \n',ASsavestr)

% PHsavestr = sprintf('%s/nt%02d_waveSet-%s_PH.mat',dirstr, nt, wp.waveSet);
% save(PHsavestr, 'ph', '-v7.3');
% fprintf('SAVED PHASE RESULTS ++++++++++ %s \n',PHsavestr)
end