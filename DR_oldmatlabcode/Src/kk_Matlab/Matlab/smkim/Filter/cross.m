% look at spectra of rawdata

load /data4/loren/The/theeeg01-2-01
raweeg1 = theeeg{1}{2}{1};
load /data4/loren/The/theeeg01-2-06
raweeg2 = theeeg{1}{2}{6};

[rpsd1 f1] = psd(raweeg1.data, 1024,  1.5024e+03);
[rpsd2 f2] = psd(raweeg2.data, 1024,  1.5024e+03);

plot(f1, rpsd1, 'b');
hold on;
plot(f2, rpsd2, 'r');
axis([0 100 0 1.5e5]);

% examin coherence
[ch f1] = cohere(raweeg1.data, raweeg2.data, 2^13, 1.5024e+03);
figure;
plot(f1, ch);
axis([0 100 0 1]);


load /data4/loren/The/thetheta01-2-01
eeg1 = thetheta{1}{2}{1};
load /data4/loren/The/thetheta01-2-06
eeg2 = thetheta{1}{2}{6};
[psd1 f1] = psd(eeg1.data, 1024,  200);
[psd2 f2] = psd(eeg2.data, 1024,  200);

%plot(f1, psd1, 'b');
%hold on;
%plot(f2, psd2, 'r');
%axis([0 100 0 1.5e5]);

% examin coherence
[fch f1] = cohere(eeg1.data, eeg2.data, 2^10, 200);
hold on
plot(f1, fch, 'r');


% phase analysis
[cr f1] = csd(raweeg1.data, raweeg2.data, 2^13, 1.5024e+03);

figure;
plot(f1(2:end-1), phase(cr(2:end-1)'));

[crtmp f1] = csd(eeg1.data, eeg2.data, 2^10, 200);
hold on
plot(f1(2:end-1), phase(crtmp(2:end-1)'), 'r');



% calculate the spectrogram
%specgram(raweeg1.data, 2^13,  1.5024e+03)
[s1 f t] = specgram(raweeg1.data, 2^13,  1.5024e+03)


[s1 f t] = specgram(eeg1.data(1:10000), 256,  200);
surf(t(1:100), f, abs(s1(1:100,:)));
