function [w, f, p] = swmts(raweeg, win)
%
%   [w, frequencies, power] = swmts(raweegstruct, win)
%
% compute sliding-window spectrogram over consecutive
% windows using multi-taper method. The output arrays are
%
%   windows, an array of the time intervals over which the 
%   spectral estimates are made. windows is an Nx2 array, 
%   where windows(:,1) contains the start times and windows(:,2) 
%   contains the end times of the windows. Note that the centers
%   of the windows can be obtained by 
%   0.5*(windows(2:end)+windows(1:end-1))
%
%   frequencies, a vector of the discrete frequency values at 
%   which spectral power is estimated. 
%
%   power, an array of the spectral power over all time 
%   time windows. Note that
%   size(power) = [ length(frequencies) length(windows) ]
%

fs= raweeg.samprate;

% discretization of frequency values for the FFT
nfft = 20*2*ceil(fs/2);

% time*bandwidth product for the tapers
nw = 5;

ndata= size(raweeg.data,1);

if nargin<2 | isempty(win)
    % size of the FFT window, in samples. The window should be
    % wide to achieve good frequency resolution, albeit at
    % the expense of temporal resolution. Here, the value
    % corresponds to a 5-second window
    windowsize = floor(5*fs); % 5 seconds

    % number of samples to step between successive windows
    % here, we are stepping in such a way that we have
    % overlapping half-windows
    stepsize = 0.5*windowsize;

    % find the earliest time which is an integer multiple of the
    % window step size, and the latest time which is an integer
    % multiple of the window step size, and count how many 
    % windows fit in between
    numwindows =  floor(length(raweeg.data)/stepsize)-1;

    % generate an array for the edges of the windows
    w(:,1) = [0:numwindows-1]'*stepsize + 1;
    w(:,2) = w(:,1) + windowsize;
else
    starttime= raweeg.starttime;
    w= floor((win-raweeg.starttime)*fs)+1;
    if find(w-ndata>10); error('inconsistency: indices messed up'); end
    w(w>ndata)= ndata;
end

numwindows= size(w,1);

% reserve an array for the output
p = zeros( nfft/2+1, numwindows );

% compute a spectral estimate for each window
%            keyboard
for i = 1:numwindows
    idxs = w(i,1):w(i,2);
    eeg= raweeg.data(idxs);
    eeg= eeg-mean(eeg);         % subtract DC component
    [p(:,i), f] = pmtm(eeg,nw,nfft,fs,'eigen');
end

% keep only the part of the spectrum that falls within 1-20 Hz
FREQUENCY_RANGE = [1 30];
keepers = find((f >= FREQUENCY_RANGE(1)) & (f <= FREQUENCY_RANGE(2)));
p = p(keepers,:);
f = f(keepers);
w= (w-1)/fs+raweeg.starttime;

%{
%imagesc(0.5*sum(w,2),f,10*log10(p));
%colorbar;
%axis xy;
%}
