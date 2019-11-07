
function [out] = dfa_eventTrigSpiking(idx, excludeperiods, varargin)

%{
run with filterframework via singlcellanal
- gathers spiking around event times (generalized from swr, lick- versions)
- example script: licktrigSUmod_20191106.m

idx: [day epoch ntrode cluster]
excludeIntervals: [start end; ...] timefilter (applies to events, spikes)
varargin required data: {<eventType>, 'spikes'}, eventType

$DR19
%}



fprintf('%d %d %d %d\n',idx)

% check for required data in varargin
reqData = {'spikes'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end

eventType = 'ca1rippleskons';
win = [0.5 0.5]; % in sec
bin = 0.001; % 1 ms for rasters
frbin= 0.01; % 10 ms for population FR plotting

if ~isempty(varargin)
    assign(varargin{:})
end

% init output
out = init_out();
out.index = idx;
day = idx(1);
ep = idx(2);
nt = idx(3);
clust = idx(4);

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
try
    eventTimes = events{day}{ep}{1}.starttime;
catch
    try
        eventTimes = events{day}{ep}.starttime;
    catch
        fprintf('no events detected for day%d ep%d\n', day, ep)
        return
    end
end

evbefore = size(eventTimes,1);
eventTimes = eventTimes(~isExcluded(eventTimes(:,1),excludeperiods),:);
evafter = size(eventTimes,1);
fprintf('%d of %d events discarded bc excluded periods in timefilter: d%d e%d\n',...
    evbefore-evafter, evbefore, day,ep)

if isempty(eventTimes)
    return
end
epStartTime = spikes{day}{ep}{nt}{clust}.timerange(1);
epEndTime = spikes{day}{ep}{nt}{clust}.timerange(2);
% Remove triggering events that are too close to the beginning or end
while eventTimes(1,1)<(epStartTime+win(1))
    eventTimes(1,:) = [];
end
while eventTimes(end,1)>(epEndTime-win(2))
    eventTimes(end,:) = [];
end

%% spiketimes
try
    spiketimes = spikes{day}{ep}{nt}{clust}.data(:,1);
catch
    spiketimes = [];
end

%% psth
time = -win(1)-0.5*bin : bin : win(2)+0.5*bin;
psth = cell2mat(arrayfun(@(r) histc(spiketimes , r + time), eventTimes(:,1), 'un', 0)')';

frtime = (-win(1)-0.5*frbin):frbin:(win(2)+0.5*frbin);
frhist = cell2mat(arrayfun(@(r) histc(spiketimes , r + frtime), eventTimes(:,1), 'un', 0)')';

instantFR = cell2mat(arrayfun(@(x) instantfr(spiketimes, x + time),eventTimes(:,1),'un',0));
%
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
out.dims = [];
out.eventTimes = [];
out.numSpikes = [];
out.time = [];
out.psth = [];

out.frtime = [];
out.frhist = [];

out.instantFR = [];
end