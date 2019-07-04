

function waveletFFT=createCMWfft(wp)
% create wavelet FFT
waveletFFT = zeros(wp.nConv2pow, wp.numfrex);
for fi=1:wp.numfrex % can use parfor
    % creating the morlet wavelet by combining the complex sine wave and the gaussian
    sine_wave = exp(2*1i*pi*wp.frex(fi).*wp.timeWin); %make a complex sine wave
    bn = wp.nWavecycles(fi)/(2*pi*wp.frex(fi)); % std of gaussian of wavelet (dependent on the freq and #cycles)
    gaus_win = exp(-wp.timeWin.^2./(2*bn^2)); % gaussian
    wavelet = sine_wave .* gaus_win;
    wavefft = fft(wavelet,wp.nConv2pow); % take the fft of the wavelet pad with zeros to the next power of 2 for speed
    % normalize wavelet to a maximum of 1 to ensure convolution units are same as data
    waveletFFT(:,fi) = (wavefft ./ max(wavefft))';
    fprintf('creating wavelet for freq %d of %d \n',fi,wp.numfrex);
end
end