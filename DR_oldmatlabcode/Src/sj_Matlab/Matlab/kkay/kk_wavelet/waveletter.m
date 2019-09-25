function [spectrogram,times,frequencies,meanspectrum,stdspec] = waveletter(data,Fs,scale,mother)

% Calls Torrence-Compo wavelet() using friendly inputs and outputs.

% Default scale is 0.25, ranges from 0 to 1.
% mother is either 'MORLET', 'PAUL','DOG'

[wave,period,scales,~] = wavelet(data,1/Fs,1,scale,5/Fs,64,'PAUL',-1);
spectrogram = (abs(wave)).^2;
frequencies = period.^(-1);

N=length(data);
meanspectrum = sum(spectrogram)/N;   % time-average over all times
spectrogram=spectrogram-repmat(meanspectrum,[N 1]);
stdspec=std(spectrogram);
spectrogram=bsxfun(@rdivide,spectrogram,stdspec);

times=(1:N-1)/Fs;

end





% MORLET  1/Fs * 3   and   3  -- ends at 250, but too smeared at lower
% DOG    1/Fs * 3 and 3 -- ends at 150, high gamma looks OK, but lower freq
            % too smeared and odd higher freq??
% MORLET 5/Fs, 64, and 6 -- fine for ripples, but bad for anything lower
    % MORLET 4/Fs, 64 and 6 -- fine for 250 Hz, but bad for anything lower
% MORLET 5/Fs, 64, and 2 -- good for high gamma?  up to 120
% MORLET 5/Fs, 64, and 12 -- up to 400+, not even good for ripples

% PAUL, 5/Fs, 64, -1  --- very good for ripples + fast gamma!
% PAUL, 4/Fs, 64, -1 -- very very good for ripples + fast gamma (decel)
% PAUL, 4/Fs, 64, 2 -- see high gamma bursts.. distorting ripples + fast g
        % max at 150
% PAUL, 4/Fs, 32, -1 -- high frequencies tinny, low frequencies blobby..
% PAUL, 4/Fs, 128, -1 -- terrible for low, ripple freq too high?
        
% DOG, 4/Fs, 64, -1
