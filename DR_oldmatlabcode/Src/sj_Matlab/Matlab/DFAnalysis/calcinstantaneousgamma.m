function out = calcinstantaneousgamma(index, excludetimes, eeg, varargin)

%Define sampling rate in Hz
smoothing_width = 0.001;
nstd = 1;
minduration = 0.003;

%Set options
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'nstd'
                nstd = varargin{option+1};
            case 'minduration'
                minduration = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

e = eeg{index(1)}{index(2)}{index(3)}.data;
times = geteegtimes(eeg{index(1)}{index(2)}{index(3)});
clear eeg

%Upsample e by a factor of 2
times = linspace(times(1),times(end),2*length(times));
e = interpft(e,length(times));
srate = 3000;

%Filter for broadband gamma
filterstring = '/home/mcarr/Src/Matlab/Filters/broadbandgammafilter3000.mat';
eval(['load ', filterstring]);
gamma = hilbert(filtfilt(gammafilter,1,e));
clear e
power = abs(gamma);
phase = angle(gamma);
%clear gamma

%Detect when gamma power exceeds nstd above the mean
% smooth the envelope:
kernel = gaussian(smoothing_width*srate, ceil(8*smoothing_width*srate));
renv = smoothvect(power, kernel);
	
%only consider periods that are greater than minduration
mindur = round(minduration * srate);

% calculate the threshold in uV units 
baseline = mean(renv);
stdev = std(renv);
thresh = baseline + nstd * stdev;

%extract the events
tmphighgamma = extractevents(renv, thresh, baseline, 0, mindur, 0)';

%convert samples to times. 
starttime = times(1) + tmphighgamma(:,1) / srate;
endtime = times(1) + tmphighgamma(:,2) / srate;

goodtimes = logical(isExcluded(times,[starttime endtime])) & ... 
    logical(~isExcluded(times,excludetimes));

%Find peaks for good times
goodphase = phase(goodtimes);
goodtime = times(goodtimes);

[maximum minimum] = peakdet(goodphase,0.1,goodtime);

%Determine the time between troughs and exclude ones that are out of range
troughs = diff([minimum(:,1); 0]);
valid = troughs>(1/140) & troughs<(1/25) & minimum(:,2)<-2;
infrequency = 1./troughs(valid);
time = minimum(valid,1);

out.gamma = infrequency;
out.time = time;
out.index = index;

end
