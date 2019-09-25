function out = calcspiketrigpower(sindex, eindex, excludetimes, eeg, spikes, varargin)
%
%

%Set Options
win = [1 1];
movingwin = [0.5 0.1];

params = {};
params.Fs = 1500;
params.fpass = [100 150];
params.trialave = 0;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendind = varargin{option+1};
            case 'win'
                win = varargin{option+1};
            case 'movingwin'
                movingwin = varargin{option+1};
            case 'fpass'
                params.fpass = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

e = eeg{eindex(1)}{eindex(2)}{eindex(3)}.data;
s = spikes{sindex(1)}{sindex(2)}{sindex(3)}{sindex(4)}.data(:,1);

% apply excludetimes
includetimes = ~isExcluded(s,excludetimes);
s = s(includetimes);

% subtract starttimes form spiketimes
s = s - eeg{eindex(1)}{eindex(2)}{eindex(3)}.starttime;

% get rid of spikes that are too close to the beginning or end
while s(1)*params.Fs < round(win(1)*params.Fs)
    s = s(2:end);
end

while (s(end)*params.Fs) > (length(e) - round(win(2)*params.Fs))
    s = s(1:end-1);
end

% calculate triggered spectrogram
[S,t,f] = mtspecgramtrigc(e,s,[0.5 0.5],[0.25 0.05],params);

out.time = t;
out.frequency = f;
out.power = S;
    
if appendind == 1
    out.index = [eindex sindex(3:4)];
end

end

