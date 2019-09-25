function phase= convertEEG2Phase(eeg, times)

peaks = findPeaks(eeg.data');
peaktimes = peaks / eeg.samprate + eeg.starttime;
phase = 2*pi*interp1(peaktimes, 1:length(peaks),times, 'linear');
valid=find(isfinite(phase));
if length(valid)/length(phase) < 0.9
    warning('Too many invalid theta phases');
end
phase=mod(phase, 2*pi);
