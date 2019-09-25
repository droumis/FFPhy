
% set the parameters for the 150 - 250 Hz bandpass filter (1500 Hz Fs)
samprate = 1500;
lpcutoff = 275;
lpfreq = 250;
hpfreq = 150;
hpcutoff = 125;

% create the f and a vectors to set the frequency response
f = [0 hpcutoff hpfreq lpfreq lpcutoff samprate/2] ./ (samprate/2);
a = [0 0 1 1 0 0];

% create the full FIR filter
fullfilt = firpm(151, f, a);

halffilt = fullfilt(1:length(fullfilt)/2);

h = spectrum.welch;
ffpsd = psd(h, fullfilt, 'Fs', samprate);
plot(ffpsd)
hold on

hfpsd = psd(h, halffilt, 'Fs', samprate);
hold on
plot(hfpsd)

