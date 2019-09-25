function out = calcspeedspectrum(index, excludetimes, eeg, pos, varargin)
%OUT = CALCEEGSPEEDSPECTRUM(INDEX,EXCLUDETIMES,EEG,POS,VARARGIN)
% computes the average spectrum for each speed bin

% Options:
% fpass: frequency to compute spectrum.
%   Default is [1 300], change this with care
% appendindex: 1 or 0, 1 to append the index.
%   Default 0.
% window: determines the time window to compute the spectrum.
%   Default is [0.5 0.5], change this with care
%

%Define options
params = {};
params.Fs = 1500;
params.fpass = [1 300];
win = [0.5 0.5];
appendindex = 0;

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'window'
                win = varargin{option+1};
            case 'appendindex'
                appendindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Initialize EEG
e = eeg{index(1)}{index(2)}{index(3)};
etimes = geteegtimes(e);
clear eeg;

% compute full spectrum
[spectrum,time,f]= mtspecgramc(e.data,[win(1) win(2)],params);
time = time+etimes(1);

% normalize full spectrum
tmp_spectrum = spectrum;
mean_spectrum = mean(spectrum);
std_spectrum = std(spectrum);
for t = 1:length(time)
    tmp_spectrum(t,:) = (spectrum(t,:)-mean_spectrum)./std_spectrum;
end
spectrum = tmp_spectrum;
clear mean_spectrum std_spectrum tmp_spectrum

% Compute mean speed for each spectral estimate
speedind_start = lookup(time-win(1)/2,pos{index(1)}{index(2)}.data(:,1));
speedind_end = lookup(time+win(1)/2,pos{index(1)}{index(2)}.data(:,1));
speed = nan(size(speedind_start));

for s = 1:length(speedind_start)
    speed(s) = mean(pos{index(1)}{index(2)}.data(speedind_start(s):speedind_end(s),8));
end

% Determine the spectrum for each speed bin
binedges = [exp(-1.75:0.1:3.25-0.25)' exp(-1.75+0.25:0.1:3.25)'];
spectrum_matrix = nan(length(f),size(binedges,1));
for i = 1:size(binedges,1)
    spectrum_matrix(:,i) = mean(spectrum(speed>binedges(i,1) & speed<binedges(i,2),:));
end

out.spectrum = spectrum_matrix;
out.bin = binedges;
out.frequency = f;

if appendindex
    out.index = index;
end

end

