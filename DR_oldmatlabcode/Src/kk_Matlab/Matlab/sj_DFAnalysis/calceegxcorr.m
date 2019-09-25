function out = calceegxcorr(index, excludetimes, eeg, varargin)
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
srate = 1500;
fpass = [100 150];
appendindex = 0;
lag = 100;
numsurrogate = 0;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'fpass'
                fpass = varargin{option+1};
            case 'lag'
                lag = varargin{option+1};
            case 'numsurrogate'
                numsurrogate = varargin{option+1};
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

% Filter for fpass
e1 = eegfilt(e1',srate,fpass(1),fpass(2));
e2 = eegfilt(e2',srate,fpass(1),fpass(2));

% Compute XCorr
[c lags] = xcorr(e1,e2,lag,'coeff');

out.xcorr = c;
out.lags = lags;

% Comput surrogate xcorr measures
if numsurrogate>0
        numpoints = length(e1);
        minskip = srate;
        maxskip = numpoints-srate;
        skip = ceil(numpoints.*rand(numsurrogate*2,1));
        skip(find(skip>maxskip))=[];
        skip(find(skip<minskip))=[];
        skip=skip(1:numsurrogate,1);
        surrogate = zeros(2*lag+1,numsurrogate);

    for s=1:numsurrogate
        surrogate_e1 = [e1(skip(s):end) e1(1:skip(s)-1)];
        surrogate(:,s) = xcorr(surrogate_e1,e2,lag);
    end

    out.surrogate = surrogate;
end

if (appendindex)
    out.index = index;
end

end