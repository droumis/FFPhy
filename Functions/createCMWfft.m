

function waveletFFT=createCMWfft(wp, nConv2pow, timeWin)
% create wavelet FFT
waveletFFT = zeros(nConv2pow, wp.numfrex);
parfor fi=1:wp.numfrex % can use parfor
    % creating the morlet wavelet by combining the complex sine wave and the gaussian
    sine_wave = exp(2*1i*pi*wp.frex(fi).*timeWin); %make a complex sine wave
    bn = wp.nWavecycles(fi)/(2*pi*wp.frex(fi)); % std of gaussian of wavelet (dependent on the freq and #cycles)
    gaus_win = exp(-timeWin.^2./(2*bn^2)); % gaussian
    wavelet = sine_wave .* gaus_win;
    wavefft = fft(wavelet,nConv2pow); % take the fft of the wavelet pad with zeros to the next power of 2 for speed
    % normalize wavelet to a maximum of 1 to ensure convolution units are same as data
    waveletFFT(:,fi) = (wavefft ./ max(wavefft))';
    fprintf('creating wavelet for freq %d of %d \n',fi,wp.numfrex);
end
end