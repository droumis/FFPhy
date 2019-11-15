
function [out] = dfa_eventTrigSpiking(idx, excludeIntervals, varargin)

%{
$version = FFPhy0.1
run with via runfilter, singleDayCellAnal
flag
- gathers spiking around event times (generalized from swr, lick- versions)
- example script: licktrigSUmod_20191106.m

idx: [day epoch ntrode cluster]
excludeIntervals: [start end; ...] timefilter (applies to events, spikes)
varargin required data: {<eventType>, 'spikes'}, eventType

$DR19
%}

% fprintf('%d %d %d %d\n',idx)

eventType = 'ca1rippleskons';
applyTFtoSpikes = 0;
win = [1 1]; % seconds. 
bin = 0.001; % seconds. rasters
% frbin= 0.01; % seconds. FR plotting
wbin = .02; % seconds. wider psth
smbins = 10; % bins. smooth across x bins (wbin x smbins = range of influence)
byDay = 1;
if ~isempty(varargin)
    assign(varargin{:})
end
fprintf('eventType %s\n', eventType)
% check for required data in varargin
reqData = {'spikes', 'cellinfo', eventType};
for s = 1:length(reqData)
    if ~any(cell2mat(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 0)))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end

day = idx(1);
if byDay
    eps = idx(4:5);
    nt = idx(2);
    clust = idx(3);
else
    eps = idx(2);
    nt = idx(3);
    clust = idx(4);
end

% init output
out = init_out();
out.index = idx;

out.cellInfo = cellinfo{day}{eps(1)}{nt}{clust};
out.area = cellinfo{day}{eps(1)}{nt}{clust}.area;
out.subarea = cellinfo{day}{eps(1)}{nt}{clust}.subarea;

%% get spikes, apply timefilter
spikeTimes = [];
numSpikesPerEp = [];
for e = 1:length(eps)
    try
        spikeTimes = [spikeTimes; spikes{day}{eps(e)}{nt}{clust}.data(:,1)];
        numSpikesPerEp = [numSpikesPerEp size(spikes{day}{eps(e)}{nt}{clust}.data,1)];
    catch
        continue
    end
end
if applyTFtoSpikes
    spikesBefore = size(spikeTimes,1);
    spikeTimes = spikeTimes(~isExcluded(spikeTimes, excludeIntervals));
    if isempty(spikeTimes)
        fprintf('spikeimes empty\n');
        return
    end
    spikesAfter = size(spikeTimes,1);
    fprintf('%d of %d spikes excluded with timefilter: day %d \n',...
        spikesBefore-spikesAfter, spikesBefore, day)
end
out.numSpikesPerEp = numSpikesPerEp;
%% get events, apply timefilter
try
    evid = find(contains(varargin(1:2:end), eventType));
    o = [1:2:length(varargin)]+1;
    events = varargin{o(evid)};
catch
    fprintf('must include event data\n');
    return
end
eventTimes = [];
numEventsPerEp = [];
evbefore = 0;
for e = 1:length(eps)
    try
        epEv = events{day}{eps(e)}{1}.starttime;
    catch
        try
            epEv = events{day}{eps(e)}.starttime;
        catch
            fprintf('no events detected for day%d ep%d\n', day, eps(e))
            continue
        end
    end
    evbefore = evbefore+size(epEv,1);
    epEvInc = epEv(~isExcluded(epEv(:,1), excludeIntervals),:);
    numEventsPerEp = [numEventsPerEp length(epEvInc)];
    eventTimes = [eventTimes; epEvInc];
end
evafter = size(eventTimes,1);
fprintf('%d of %d events excluded with timefilter: d%d \n',...
    evbefore-evafter, evbefore, day)

if isempty(eventTimes)
    fprintf('eventTimes is empty\n');
    return
end

%% Remove events that are too close to the beginning or end
epStartTime = [];
epEndTime = [];
for e = 1:length(eps)
    epStartTime = [epStartTime spikes{day}{eps(e)}{nt}{clust}.timerange(1)];
    epEndTime = [epEndTime spikes{day}{eps(e)}{nt}{clust}.timerange(2)];
end
epStartTime = min(epStartTime);
epEndTime = max(epEndTime);

while eventTimes(1,1)<(epStartTime+win(1))
    eventTimes(1,:) = [];
end
while eventTimes(end,1)>(epEndTime-win(2))
    eventTimes(end,:) = [];
end

%% psth
time = -win(1)-0.5*bin : bin : win(2)+0.5*bin;
wtime = -win(1)-0.5*wbin : wbin : win(2)+0.5*wbin;
% frtime = -win(1)-0.5*frbin: frbin: win(2)+0.5*frbin;

% spike psth
psth = cell2mat(arrayfun(@(r) histc(spikeTimes , r + time), eventTimes(:,1), 'un', 0)')';
psthWZ = zscore(cell2mat(arrayfun(@(r) histcounts(spikeTimes , r + wtime), eventTimes(:,1), 'un', 0)),[],2);
psthWZM = nanmean(psthWZ);
psthWZMS = smoothdata(psthWZM, 'loess', smbins);

% instantaneous firing rate 
psfr = cell2mat(arrayfun(@(x) instantfr(spikeTimes, x + time),eventTimes(:,1),'un',0));
psfrWZ = zscore(cell2mat(arrayfun(@(x) instantfr(spikeTimes, x + wtime), eventTimes(:,1),'un',0)),[],2);
psfrWZM = nanmean(psfrWZ);
psfrWZMS = smoothdata(psfrWZM, 'loess', smbins);

% binned firing rate
% frhist = cell2mat(arrayfun(@(r) histc(spiketimes , r + frtime), eventTimes(:,1), 'un', 0)')';

%%
out.numEventsPerEp = numEventsPerEp;
out.numSpikesPerEp = numSpikesPerEp;
out.dims = {'event', 'time'};
out.eventTimes = eventTimes;
out.numSpikes = length(spikeTimes);

out.time = time;
out.wtime = wtime;

out.psth = psth;
out.psthWZ = psthWZ;
out.psthWZM = psthWZM;
out.psthWZMS = psthWZMS;

out.psfr = psfr;
out.psfrWZ = psfrWZ;
out.psfrWZM = psfrWZM;
out.psfrWZMS = psfrWZMS;

end

function out = init_out()
out.index = [];
out.cellInfo = []; 
out.area = '';
out.subarea = '';
out.numEventsPerEp = [];
out.numSpikesPerEp = [];
out.dims = [];
out.eventTimes = [];
out.numSpikes = [];

out.time = [];
out.wtime = [];

out.psth = [];
out.psthWZ = [];
out.psthWZM = [];
out.psthWZMS = [];

out.psfr = [];
out.psfrWZ = [];
out.psfrWZM = [];
out.psfrWZMS = [];
end