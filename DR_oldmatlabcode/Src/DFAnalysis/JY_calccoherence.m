
function out = calccoherence(index, excludetimes, data, eeg, varargin)
% function out = calccoherence_test(index, excludetimes, eeg, varargin)
%
%  Plots the coherence for an eeg tetrode pair. If you use a time filter,
%  excluded times are removed and the includedtimes are averaged together.
%
%   out is a structure with the following fields
%       coherence-- This is the coherence for the tetrode pairs
%       frequency-- Frequency vector
%       index-- Only if appendindex is set to 1 (default)


% parse the options
params = [];
params.Fs = 1500;
params.fpass = [1 250];
params.err = [1 0.05];
appendindex = 1;
window =[10 0.5];
params.tapers=[3 5];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'Fs'
                params.Fs = varargin{option+1};
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'err'
                params.err = varargin{option+1};
	    otherwise
                error(['Option ',varargin{option},' unknown.']);
        end   
    else
        error('Options must be strings, followed by the variable');
    end
end


% assign a temporary variable for eeg
e1 = eeg{index(1)}{index(2)}{index(3)};
e2 = eeg{index(1)}{index(2)}{index(4)};
e1start = e1.starttime;

e1times = geteegtimes(e1);
e2times = geteegtimes(e2);

if length(e1times)>length(e2times)
    temp = lookup(e2times,e1times);
    e1 = e1.data(temp);
    e2 = e2.data;
elseif length(e2times)>length(e1times)
    temp = lookup(e1times,e2times);
    e1 = e1.data;
    e2 = e2.data(temp);
elseif length(e1times)==length(e2times)
    e1 = e1.data;
    e2 = e2.data;
end

% compute full coherence
[C,phi,S12,S1,S2,t,f]=cohgramc(e1,e2,window, params);

%apply excludetimes
goodtimes = ~isExcluded(t+e1start, excludetimes);
tempcohere = C;
tempcohere= tempcohere(goodtimes,:);
coherence = tempcohere;%mean(tempcohere);

p = data{index(1)}{index(2)}.Pos.correcteddata;
time = lookup(t+e1start, p(:,1));
speed = p(time,5);

out.speed = speed(goodtimes);
out.coherence = coherence;
out.frequency = f;

if (appendindex)
    out.index = index;
end

end
