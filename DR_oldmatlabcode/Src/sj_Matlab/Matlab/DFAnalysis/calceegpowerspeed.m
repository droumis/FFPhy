function out = calceegpowerspeed(index, excludetimes, eeg, pos,linpos, varargin)
%OUT = CALCEEGPOWERSPEED(INDEX,EXCLUDETIMES,EEG,POS,VARARGIN)
%
% calceegpowerspeed computes a windowed spectrogram and for each time
% computes the speed, theta frequency, and theta:delta power ratio.
%
% out.power     theta:delta ratio
% out.frequency dominant frequency in the theta range
% out.speed     estimate at speed for each time computed using a weighted
%               average of speeds surrounding the center using a 15-sample
%               wide hann window.

%Define options
params = {};
params.Fs = 1500;
params.fpass = [1 300];
win = [0.5 0.5];
appendindex = 1;

delta_frequency = [2 4];
theta_frequency = [6 10];
lowgamma_frequency = [20 55];
highgamma_frequency = [65 140];

cellfilter = [];

if ~isempty(varargin)
    if isstruct(varargin{1}(1))
        baseline = varargin{1}(1).mean;
        stdev = varargin{1}(1).std;
        varargin(1) = [];
    end
end

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'window'
                win = varargin{option+1};
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'delta'
                delta_frequency = varargin{option+1};
            case 'theta'
                theta_frequency = varargin{option+1};
            case 'lowgamma'
                lowgamma_frequency = varargin{option+1};
            case 'highgamma'
                highgamma_frequency = varargin{option+1};
            case 'cellfilter'
                cellfilter = varargin{option+1};
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

if time > 1000
    cutoff = 100;
else
    cutoff = 50;
end

% normalize full spectrum
time = time(cutoff:end-cutoff);
spectrum = spectrum(cutoff:end-cutoff,:);
tmp_spectrum = spectrum;
mean_spectrum = mean(spectrum);
std_spectrum = std(spectrum);

for t = 1:length(time)
    if exist('baseline','var')
        tmp_spectrum(t,:) = (spectrum(t,:) - baseline)./stdev;
    else
        tmp_spectrum(t,:) = (spectrum(t,:)-mean_spectrum)./std_spectrum;
    end
end
spectrum = tmp_spectrum;

% Identify frequencies of theta, delta, low gamma, and high gamma
delta = lookup([delta_frequency(1) delta_frequency(2)],f);
theta = lookup([theta_frequency(1) theta_frequency(2)],f);
lowgamma = lookup([lowgamma_frequency(1) lowgamma_frequency(2)],f);
highgamma = lookup([highgamma_frequency(1) highgamma_frequency(2)],f);

% Compute the theta:delta ratio and theta power
theta_power = sum(spectrum(:,theta(1):theta(2)),2);
delta_power = sum(spectrum(:,delta(1):delta(2)),2);
theta_delta_ratio = theta_power./delta_power;

%Compute low and high gamma power
lowgamma_power = sum(spectrum(:,lowgamma(1):lowgamma(2)),2);
highgamma_power = sum(spectrum(:,highgamma(1):highgamma(2)),2);

% Compute mean speed at each time
speedind_start = lookup(time-win(1)/2,pos{index(1)}{index(2)}.data(:,1));
speedind_end = lookup(time+win(1)/2,pos{index(1)}{index(2)}.data(:,1));
speed = nan(size(speedind_start));
linearind_start = lookup(time-win(1)/2,linpos{index(1)}{index(2)}.statematrix.time);
linearind_end = lookup(time+win(1)/2,linpos{index(1)}{index(2)}.statematrix.time);
linearspeed = nan(size(linearind_start));

for s = 1:length(speedind_start)
    speed(s) = mean(pos{index(1)}{index(2)}.data(speedind_start(s):speedind_end(s),8));
    linearspeed(s) = mean(abs(linpos{index(1)}{index(2)}.statematrix.linearVelocity(linearind_start(s):linearind_end(s),1)));
end

out.theta_delta_ratio = theta_delta_ratio;
out.theta_power = theta_power;
out.delta_power = delta_power;
out.speed = speed;
out.linearspeed = linearspeed;
out.lowgamma_power = lowgamma_power;
out.highgamma_power = highgamma_power;
out.mean = mean_spectrum;
out.std = std_spectrum;
out.time = time;

if appendindex
    out.index = index;
end

end