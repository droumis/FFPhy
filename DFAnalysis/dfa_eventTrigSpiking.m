
function [out] = dfa_eventTrigSpiking(idx, excludeIntervals, varargin)

%{
run with filterframework via singlcellanal or singleDayCellAnal via byDay
flag
- gathers spiking around event times (generalized from swr, lick- versions)
- example script: licktrigSUmod_20191106.m

idx: [day epoch ntrode cluster]
excludeIntervals: [start end; ...] timefilter (applies to events, spikes)
varargin required data: {<eventType>, 'spikes'}, eventType

$DR19
%}

% fprintf('%d %d %d %d\n',idx)

% check for required data in varargin
reqData = {'spikes'};
for s = 1:length(reqData)
    if ~any(cell2mat(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 0)))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end

eventType = 'ca1rippleskons';
win = [1 1]; % in sec
bin = 0.001; % 1 ms for rasters
frbin= 0.01; % 10 ms for population FR plotting
byDay = 1;
if ~isempty(varargin)
    assign(varargin{:})
end

% init output
out = init_out();
out.index = idx;
day = idx(1);
if byDay
    ep = idx(4:5);
    nt = idx(2);
    clust = idx(3);
else
    ep = idx(2);
    nt = idx(3);
    clust = idx(4);
end
%% get events from ~excludeperiods
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
for e = 1:length(ep)
    try
        epEv = events{day}{ep(e)}{1}.starttime;
    catch
        try
            epEv = events{day}{ep(e)}.starttime;
        catch
            fprintf('no events detected for day%d ep%d\n', day, ep(e))
            continue
        end
    end
    epEvInc = epEv(~isExcluded(epEv(:,1), excludeIntervals),:);
    eventTimes = [eventTimes; epEvInc];
    numEventsPerEp = [numEventsPerEp length(epEvInc)];
end
% evbefore = size(eventTimes,1);
% evafter = size(eventTimes,1);
% fprintf('%d of %d events discarded bc excluded periods in timefilter: d%d \n',...
%     evbefore-evafter, evbefore, day)

if isempty(eventTimes)
    fprintf('eventTimes is empty\n');
    return
end
epStartTime = [];
epEndTime = [];
for e = 1:length(ep)
    epStartTime = [epStartTime spikes{day}{ep(e)}{nt}{clust}.timerange(1)];
    epEndTime = [epEndTime spikes{day}{ep(e)}{nt}{clust}.timerange(2)];
end
epStartTime = min(epStartTime);
epEndTime = max(epEndTime);
% Remove triggering events that are too close to the beginning or end
while eventTimes(1,1)<(epStartTime+win(1))
    eventTimes(1,:) = [];
end
while eventTimes(end,1)>(epEndTime-win(2))
    eventTimes(end,:) = [];
end

%% spiketimes
spiketimes = [];
numSpikesPerEp = [];
for e = 1:length(ep)
    try
        spiketimes = [spiketimes; spikes{day}{ep(e)}{nt}{clust}.data(:,1)];
        numSpikesPerEp = [numSpikesPerEp size(spikes{day}{ep(e)}{nt}{clust}.data,1)];
    catch
        continue
    end
end
%% psth
time = -win(1)-0.5*bin : bin : win(2)+0.5*bin;
psth = cell2mat(arrayfun(@(r) histc(spiketimes , r + time), eventTimes(:,1), 'un', 0)')';

frtime = (-win(1)-0.5*frbin):frbin:(win(2)+0.5*frbin);
frhist = cell2mat(arrayfun(@(r) histc(spiketimes , r + frtime), eventTimes(:,1), 'un', 0)')';

instantFR = cell2mat(arrayfun(@(x) instantfr(spiketimes, x + time),eventTimes(:,1),'un',0));
%
out.numEventsPerEp = numEventsPerEp;
out.numSpikesPerEp = numSpikesPerEp;
out.dims = {'event', 'time'};
out.eventTimes = eventTimes;
out.numSpikes = length(spiketimes);
out.time = time;
out.psth = psth;

out.frtime = frtime;
out.frhist = frhist;

out.instantFR = instantFR;

% out.posteventmatrix = posteventmatrix;
% out.eventduration = eventduration;
% out.nospikes = nospikes;
% out.noevents = noevents;
end

function out = init_out()
out.index = [];
out.numEventsPerEp = [];
out.numSpikesPerEp = [];
out.dims = [];
out.eventTimes = [];
out.numSpikes = [];
out.time = [];
out.psth = [];

out.frtime = [];
out.frhist = [];

out.instantFR = [];
end