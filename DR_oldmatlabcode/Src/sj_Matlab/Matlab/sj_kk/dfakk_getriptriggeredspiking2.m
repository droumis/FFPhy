function [out] = dfakk_getriptriggeredspiking2(index, excludeperiods, ripple, ripples, spikes, tetinfo, varargin)
% [out] = getriptriggeredspiking(index, excludeperiods, ripples, spikes, varargin)
%
%   This function creates a ripple triggered spiking average. The ripples
%   are defined by the first group of cells and the spikes are from the
%   second group of cells.
%
%   index [day epoch tetrode cell tetrode cell]
%
%   options are
%	'minthresh',
%		     specifies the minimum threshold of a valid ripple event
%   'window', 1x2 vector specifies the window before and after each included ripple.
%                   Default is 100 mseconds before and 15 seconds after
%                   ripple start time.
%   'tetfilter' -- use this to filter for tetrodes for which to find
%                    consensus ripples

%   out = out.out   An R x C sized matrix where each entry is the number of
%                   spikes cell C fired during ripple R
%         out.index [D E T C], gives the identity of the cells for each
%                   column in out.out
%         out.times [starttime endtime], givest the starttime and endtime
%                   of the ripples for each row in out.out

% assign the options
minthresh = 0;
window = [0.1 0.1];
binsize = 0.01;

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'minthresh'
            minthresh = varargin{option+1};
        case 'window'
            window = varargin{option+1};
        case 'tetfilter'
            tetfilter = varargin{(option+1)};   % usual tetfilter string : '(isequal($area, ''CA1''))'
        case 'binsize'
            binsize = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end



% First receive valid ripple periods for the epoch's day.
% Second for each include period, find all the participating tetrodes' time of peak ripple power. 

day = index(1,1);
epoch = index(1,2);

if ~isempty(tetfilter)
    tetlist =  evaluatefilter(tetinfo{day}{epoch},tetfilter); 
else
    % default is ALL tetrodes w/ extracted ripples
    tetlist = find(cellfun(@isempty,ripples{day}{epoch})==0); 
end

r = ripples{day}{epoch}{tetlist(1)};
times = r.timerange(1):0.001:r.timerange(end);         % 1 ms step vector
%obtain includetimes
includetimes = ~isExcluded(times, excludeperiods);     % list of ones and zeros sampled every millisecond, ones = included, zeros = excluded

if sum(includetimes)==0
    disp(sprintf('no includetimes, day %d epoch %d',day,epoch));
    disp(sprintf('%d tetrodes this epoch',length(tetlist)));
    out.cellindices = [];
    out.time = [];
    out.c1vsc2 = [];
    out.nospikes = [];
    out.noripples = [];
    return
end


%reconstitute valid ripple into rippleperiod
ripstart = times([0 diff(includetimes')]==1);
ripend = times(diff(includetimes')==-1);
riptimes = [ripstart' ripend'];

rippleperiods = struct;

% iterate through each valid rippleperiod to find participating tetrodes

for p=1:size(riptimes,1)
    
    rippleperiods(p).startrip = riptimes(p,1);     % time in 1 ms vector
    rippleperiods(p).endrip = riptimes(p,2);
    rippleperiods(p).tetrodes = [];
    rippleperiods(p).peakpowers = [];
    rippleperiods(p).peaktimes = [];
    rippleperiods(p).starttimes = [];
    rippleperiods(p).waveforms = [];
    flag = 0;
    
    for t=1:length(tetlist)
        ind = find((ripples{day}{epoch}{tetlist(t)}.midtime > riptimes(p,1)) & ... 
                   (ripples{day}{epoch}{tetlist(t)}.midtime < riptimes(p,2)))' ;
        if length(ind) > 1
            disp(sprintf('tet %d -- multi ripple! (%d extracted rips in includeperiod)',tetlist(t),length(ind)))
        end
        for jj=ind       % iteration here since a tetrode may have two detected periods in the includeperiod
            rippleperiods(p).tetrodes = [rippleperiods(p).tetrodes tetlist(t)];
            rippleperiods(p).peakpowers = [rippleperiods(p).peakpowers ripples{day}{epoch}{tetlist(t)}.peak(jj)];
            rippleperiods(p).peaktimes = [rippleperiods(p).peaktimes ripples{day}{epoch}{tetlist(t)}.midtime(jj)];
            rippleperiods(p).starttimes = [rippleperiods(p).starttimes ripples{day}{epoch}{tetlist(t)}.starttime(jj)];
                if flag == 0     % this block ensures waveforms of the same vector length across all tetrodes
                    ripeegtimes = geteegtimes(ripple{day}{epoch}{tetlist(t)});
                    waveformstart = lookup(riptimes(p,1),ripeegtimes);
                    waveformend = lookup(riptimes(p,2),ripeegtimes);
                    waveform_nosamples = waveformend-waveformstart+1;
                    flag = 1;
                end
            rippleperiods(p).waveforms = [rippleperiods(p).waveforms ; ...
                                          ripple{day}{epoch}{tetlist(t)}.data(waveformstart:(waveformstart+waveform_nosamples),3)'];
        end
        
    end
    
    rippleperiods(p).averagepeaktime = mean(rippleperiods(p).peaktimes);
    rippleperiods(p).maxpeakpower = max(rippleperiods(p).peakpowers);
    
end


%% mini-study, the variance of ripple start times and peak times

if 0
start_devs1=[]; start_devs2=[];
peak_devs=[]; start_vars=[];
peak_vars=[];

for p=1:length(rippleperiods) 
    start_vars = [start_vars var(rippleperiods(p).starttimes)];
    peak_devs = [peak_devs var(rippleperiods(p).peaktimes)];
end

figure
hist(start_vars,100);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[1 0 0],'EdgeColor','w')
[N,X] = hist(start_vars,100);
hold on
alpha(0.5)
hist(peak_devs,X)
%set(gca,'yscale','log')
end

%%

% select which event parameter to use, collect them

riptimes = [];

for p = 1:length(rippleperiods)
    riptimes = [riptimes ; rippleperiods(p).startrip];
    %riptimes = [ripstart ; rippleperiods(p).averagepeaktime];
end




%Calculate the xcorrelation between spikes and ripple starttimes.

for i = 1:size(index,1)     % iterate over cells
    spiketimes = spikes{index(i,1)}{index(i,2)}{index(i,3)}{index(i,4)}.data(:,1);
    if ~isempty(riptimes)
        xcorrstruct = spikexcorr(spiketimes, riptimes, binsize, window(2));
        c1vsc2(i,:) = xcorrstruct.c1vsc2;       % unnormalized cross-correlogram
        nospikes(i,:) = xcorrstruct.nspikes1;
        noripples(i,:) = xcorrstruct.nspikes2;
        cellindices(i,:) = index(i,:);
    end
end

if exist('c1vsc2')
    
out.cellindices = cellindices;
out.time = xcorrstruct.time;
out.c1vsc2 = c1vsc2;
out.nospikes = nospikes;
out.noripples = noripples;

end
end