
function [out] = dfa_eventTrigSpiking(idx, timeFilter, varargin)
% [out] = dfa_eventTrigSpiking(idx, timeFilter, varargin)
% Gather Spikes around event times
%               
%                   Redolent Rat
%           ,,==.
%          //    `
%         ||      ,--~~~~-._ _(\--,_
%          \\._,-~   \      '    *  `o
%           `---~\( _/,___( /_/`---~~
%                 ``==-    `==-,
% Iterator:
% - singleDayCellAnal
% 
% args:
% idx: [day epoch ntrode cluster]
% timeFilter: intervals to exclude. [start end; ...] applies to events, spikes
% 
% varargs:
% - data: (i.e. 'spikes', spikes)
% - win:
% - eventType:
% - applyTFtoSpikes:
% - win: seconds. 
% - bin: seconds.
% - wbin: seconds. wider psth
% - smbins: bins. smooth across bins
%{
Notes:
- barn:rat:wheelbarrow
- gathers spiking around event times (generalized from swr, lick- versions)
- probably outdated example notebook: licktrigSUmod_20191106.m

FFPhy V0.1
@DKR
%}

eventType = 'ca1rippleskons';
applyTFtoSpikes = 0;
win = [1.5 1.5]; % seconds. 
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
    spikeTimes = spikeTimes(~isExcluded(spikeTimes, timeFilter));
    if isempty(spikeTimes)
        fprintf('spikeimes empty\n');
        return
    end
    spikesAfter = size(spikeTimes,1);
    fprintf('spikes excluded by timeFilter: %d of %d spikes\n',...
        spikesBefore-spikesAfter, spikesBefore)
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
    epEvInc = epEv(~isExcluded(epEv(:,1), timeFilter),:);
    numEventsPerEp = [numEventsPerEp length(epEvInc)];
    eventTimes = [eventTimes; epEvInc];
end
evafter = size(eventTimes,1);
fprintf('events excluded by timefilter: %d of %d\n',...
    evbefore-evafter, evbefore)

if isempty(eventTimes)
    fprintf('eventTimes is empty\n');
    return
end

% Remove events that are too close to the beginning or end
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

%% stack the event trig spikes
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