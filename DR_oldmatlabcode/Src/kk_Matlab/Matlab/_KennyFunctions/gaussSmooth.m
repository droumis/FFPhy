function w = gaussSmooth(Fs,sigma,L)
% function w = gaussSmooth(Fs,sigma,L)

% Fs - sampling rate (Hz)
% sigma - standard deviation of gaussian
% L - number of standard deviations (single side)

n = -round(sigma*Fs*L):round(sigma*Fs*L); % t = n/Fs

w = exp(-0.5 * (n / Fs/ sigma).^2);

w = w./sum(w);
