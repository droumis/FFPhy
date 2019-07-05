function saveAS(as,nt,ntN,wp, animal, lfptype)

pathdef = animaldef(lower('Demetris'));
dirstr = sprintf('%s/analyticSignal/%s', pathdef{2}, animal);
if ~isdir(dirstr)
    mkdir(dirstr);
end

analyticsignal.as = as;
analyticsignal.waveparams = wp;
analyticsignal.ntrodeN = ntN;
analyticsignal.ntrode = nt;
analyticsignal.lfptype = lfptype;

ASsavestr = sprintf('%s/AnalyticSignal_nt%02d_waveSet-%s_%s.mat', dirstr, ...
    analyticsignal.ntrode, analyticsignal.waveparams.waveSet, lfptype);
save(ASsavestr, 'analyticsignal', '-v7.3');
fprintf('SAVED ANALYTIC SIGNAL RESULTS ++++++++++ %s \n',ASsavestr)

% PHsavestr = sprintf('%s/nt%02d_waveSet-%s_PH.mat',dirstr, nt, wp.waveSet);
% save(PHsavestr, 'ph', '-v7.3');
% fprintf('SAVED PHASE RESULTS ++++++++++ %s \n',PHsavestr)
end