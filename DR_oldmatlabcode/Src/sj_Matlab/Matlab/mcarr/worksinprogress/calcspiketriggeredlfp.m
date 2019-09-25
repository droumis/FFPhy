function out = calcspiketriggeredlfp(index, excludeperiods, eeg, ripples, spikes, cellinfo, varargin)
% function out = calcspiketriggeredlfp(index, excludeperiods, eeg, ripples, spikes, cellinfo varargin)
%
% For each cell this function pulls out the relevant lfp signal at the time
% of each spike +/- bin and calculates the approximate phase in degrees.

% set options
bin = 0.01; %default time bin is 10ms on either side of a spike
minthresh = 3; %default ripple threshold is 3
cellfilter = [];

out.phase = {}; out.ripstd = {};

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'bin'
                bin = varargin{option+1};
             case 'minthresh'
                minthresh = varargin{option+1};
            case 'cellfilter'
                cellfilter = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% assign variables
e = eeg{index(1)}{index(2)}{index(3)}; clear eeg;
s = spikes{index(1)}{index(2)}{index(3)}; clear spikes;
times = geteegtimes(e);
timestep = median(diff(times));

%Find valid ripples
if isempty(cellfilter)
    [riptimes ripstd] = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        '(isequal($area, ''CA1''))', 'excludeperiods', excludeperiods,'minstd',minthresh);
else
    [riptimes ripstd] = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        cellfilter,'excludeperiods', excludeperiods,'minstd',minthresh);
end
if ~isempty(riptimes)
    %Exclude ripples occuring close together for preceding modulationg
    valid_ripples = [1000; riptimes(2:end,2)-riptimes(1:end-1,1)];
    valid_ripples = valid_ripples > 1;
    riptimes = riptimes(valid_ripples,:); clear valid_ripples;
end

%Identify valid cells on this tetrode
cellindex = evaluatefilter(cellinfo{index(1,1)}{index(1,2)}{index(1,3)},'$numspikes > 100 & $meanrate < 7');

%Pre-allocate output
out.ripstd = cell(length(cellindex),1);
out.lfp = cell(length(cellindex),1);
zero_time = round(bin/timestep)+1;
%Go through each cell and pull out the spike triggered lfp +/- bin
if ~isempty(riptimes) && ~isempty(cellindex)
    for c = cellindex'
        spiketimes = s{c}.data(:,1);
        
        %Identify spikes that occur near ripples
        spikebins = periodAssign(spiketimes,riptimes(:,[1 2]));
        spike_ripstd = ripstd(spikebins(spikebins>0));
        spiketimes = spiketimes(spikebins > 0); clear spikebins
        
        if ~isempty(spiketimes)
            %Exclude spikes too near the edge
            while spiketimes(1) < times(1)-bin
                spiketimes(1) = [];
                spike_ripstd(1) = [];
            end
            while spiketimes(end) > times(end)-bin
                spiketimes(end) = [];
                spike_ripstd(end) = [];
            end
            xcorr_rip = zeros(length(spiketimes),length(-round(bin/timestep):1:round(bin/timestep)));
            out.ripstd{c==cellindex} = spike_ripstd; clear spike_ripstd
            if ~isempty(spiketimes)
                spikeind = lookup(spiketimes,times);
                phase = nan(size(spikeind));
                for i = 1:length(spikeind)
                    %Pull out windows around each spike
                    tmp = e.data(spikeind(i)-round(bin/timestep):spikeind(i)+round(bin/timestep),1);
                    tmp_time = times(spikeind(i)-round(bin/timestep):1:spikeind(i)+round(bin/timestep));
                    xcorr_rip(i,:) = tmp;
%                     %Determine local maxima
%                     [maxtab mintab] = peakdet(tmp,0.5);
%                     previous_peak = find(maxtab(:,1)<zero_time,1,'last');
%                     next_peak = find(maxtab(:,1)>= zero_time,1,'first');
%                     
%                     %Create linear spaced phase between these two points
%                     tmp_time = linspace(tmp_time(maxtab(previous_peak,1)),tmp_time(maxtab(next_peak,1)),360);
%                     
%                     phase(i) = lookup(spiketimes(i),tmp_time);
%                     %clear mintab previous_peak next_peak tmp
                end
                %out.phase{c==cellindex} = phase; clear tmp
                out.lfp{c==cellindex} = xcorr_rip;
            end
        end
    end
end

end
