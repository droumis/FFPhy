function out = calcinstantaneousfrequency(index, excludetimes, eeg, ripples,cellinfo,varargin)

%Define sampling rate in Hz
minthresh = 3;
cellfilter = [];

%Set options
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'cellfilter'
                cellfilter = varargin{option+1};
            case 'minthresh'
                mintrhesh = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Find valid ripples
if isempty(cellfilter)
    riptimes = getripples([index(1) index(2)], ripples, cellinfo, 'cellfilter', ...
        '(isequal($area, ''CA1''))','excludeperiods', excludetimes,'minstd',minthresh);
else
    riptimes = getripples([index(1) index(2)], ripples, cellinfo, 'cellfilter', ...
        cellfilter,'excludeperiods', excludetimes,'minstd',minthresh);
end

%Exclude ripples occuring close together for preceding modulationg
valid_ripples = [1000; riptimes(2:end,2)-riptimes(1:end-1,1)];
valid_ripples = valid_ripples > 1;
riptimes = riptimes(valid_ripples,:);

if ~isempty(riptimes)
    %Filter for beta-gamma
    e = eeg{index(1)}{index(2)}{index(3)}.data;
    times = geteegtimes(eeg{index(1)}{index(2)}{index(3)});
    clear eeg valid_ripples ripples cellinfo mintrhesh excludetimes

    filterstring = '/home/mcarr/Src/Matlab/Filters/betagammafilter.mat';
    eval(['load ', filterstring]);
    gamma = hilbert(filtfilt(betagammafilter,1,e));
    phase = angle(gamma);
    clear e gamma filterstring
    
    % Go through each valid ripple and determine the histogram of frequency
    rip_start = lookup(riptimes(:,1),times);
    rip_end = lookup(riptimes(:,2),times);
    infrequency = []; ripple = [];
    for r = 1:size(rip_start,1)
        [maximum minimum] = peakdet(phase(rip_start(r):rip_end(r)),0.5,times(rip_start(r):rip_end(r)));
        if ~isempty(minimum)
            troughs = diff([minimum(:,1); times(rip_end(r))]);
            valid = troughs>=(1/50) & troughs<=(1/10) & minimum(:,2)<-2;
            if sum(valid)>3
                infrequency = [infrequency; 1./troughs(valid)];
                ripple = [ripple; r*ones(sum(valid),1)];
            end
        end
    end
	out.freq = infrequency;
    out.rip = ripple;
else
    out = [];
end

end
