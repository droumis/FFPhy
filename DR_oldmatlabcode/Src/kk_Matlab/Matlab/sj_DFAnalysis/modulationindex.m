function out = modulationindex(index, excludetimes, eeg, varargin)

%Define default options
nbin = 10;
frequency = 20:2:140;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'nbin'
                nbin = varargin{option+1};
            case 'frequency'
                frequency = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

mi = zeros(length(frequency),1);
modulationindex = zeros(length(frequency),2*nbin-1);

e = -1*eeg{index(1)}{index(2)}{index(3)}.data;
times = geteegtimes(eeg{index(1)}{index(2)}{index(3)});

clear eeg

%Compute theta phase at each time for the phase offset
load '/home/mcarr/Src/Matlab/Filters/thetafilter.mat'
theta = filtfilt(thetafilter,1,e);
phase = angle(hilbert(theta));
    
%Find valid theta intervals
goodtimes = getvalideegtimes(phase,times);
goodtimes = logical(goodtimes) & ~isExcluded(times,excludetimes);
phase = phase(logical(goodtimes));
theta = theta(logical(goodtimes));
center = (-pi+pi/(2*nbin)):pi/nbin:(pi-pi/(2*nbin));
phase = lookup(phase,center(2:end));
temptheta = accumarray(phase,theta,[],@(x) mean(x));

for g = frequency(1:end)
    r = find(frequency==g);
    %Compute analytic amplitude time series
    filterstring = sprintf('/home/mcarr/Src/Matlab/Filters/cfcampfilt%s.mat',num2str(g));
    eval(['load ', filterstring]);
    amplitude = abs(hilbert(filtfilt(cfcampfilt,1,e)));
    amplitude = amplitude(logical(goodtimes));
    %Compute histogram of gamma amplitude over theta phase for all valid times
    temp = accumarray(phase, amplitude, [2*nbin-1 1], @(x) mean(x));
    modulationindex(r,:) = temp/sum(temp);
    mi(r) = (log(2*nbin-1) + sum(modulationindex(r,:).*log(modulationindex(r,:))))/log(2*nbin-1);
end

out.modulationindex = modulationindex;
out.mi = mi;
out.theta = temptheta;
out.index = index;
end