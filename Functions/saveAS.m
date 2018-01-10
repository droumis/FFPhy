function saveAS(as, ph, dirstr, nt, wp)

ASsavestr = sprintf('%s/nt%02d_waveSet-%s_AS.mat',dirstr, nt, wp.waveSet);
save(ASsavestr, 'as', '-v7.3');
disp(sprintf('SAVED ANALYTIC SIGNAL RESULTS ++++++++++ %s',ASsavestr))
PHsavestr = sprintf('%s/nt%02d_waveSet-%s_PH.mat',dirstr, nt, wp.waveSet);
save(PHsavestr, 'ph', '-v7.3');
disp(sprintf('SAVED PHASE RESULTS ++++++++++ %s',PHsavestr))
end