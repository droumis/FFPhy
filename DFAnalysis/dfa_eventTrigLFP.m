function out = dfa_eventTrigLFP(idx, timeFilter, varargin)
% out = dfa_eventTrigLFP(idx, excludeIntervals, varargin)
% Gather LFP around event times
%
%                 Bellicose Bear
%             .--.              .--.
%            : (\ ". _......_ ." /) :
%             '.    `        `    .'
%              /'   _        _   `\
%             /     0}      {0     \
%            |       /      \       |
%            |     /'        `\     |
%             \   | .  .==.  . |   /
%              '._ \.' \__/ './ _.'
%              /  ``'._-''-_.'``  \
%Iterator:
% - multitetrodeanal
% timeFilter: reconstructs events (CHANGE THIS)
% 
% args:
% - idx [day epoch]
% - excludeIntervals
% 
% varargs:
% - data: (i.e. 'eeg', eeg)
% - win:
% - eventType:
% - LFPTypes: 
%{

Notes:
- forest:bear:cactus:mushroom:leaf

FFPhy V0.1
@DR
%}


eventType = 'lick';
LFPtypes = {'eeg', 'ripple'};
win = [0.5 0.5];
if ~isempty(varargin)
    assign(varargin{:});
end

fprintf('eventType %s\n', eventType)
% check for required data in varargin
reqData = {LFPtypes{:}, eventType};
for s = 1:length(reqData)
    if ~any(cell2mat(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 0)))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end

day = idx(1,1);
ep = idx(1,2);
nts = idx(:,3);
win = abs(win);
% init output
out = init_out();
% out.ntinfo = tetinfo{day}{ep}{nt};

%% get LFP

lfpTypeIdx = find(contains(varargin(1:2:end), LFPtypes{1}));
o = [1:2:length(varargin)]+1;
e = varargin{o(lfpTypeIdx)};
num_samp = length(e{day}{ep}{nts(end)}.data);
samprate = e{day}{ep}{nts(end)}.samprate;
epStartTime = double(e{day}{ep}{nts(end)}.starttime);
epEndTime = (num_samp/samprate) + epStartTime;
w = win*samprate;

%% get event times, apply timefilter
evid = find(contains(varargin(1:2:end), eventType));
events = varargin{o(evid)};
eventTimes = [];
try
    eventTimes = events{day}{ep}{1}.starttime;
catch
    try
        eventTimes = events{day}{ep}.starttime;
    catch
        fprintf('no events detected for day%d ep%d\n', day,ep)
        return
    end
end

evbefore = size(eventTimes,1);
eventTimes = eventTimes(~isExcluded(eventTimes(:,1),timeFilter),:);
evafter = size(eventTimes,1);
fprintf('events excluded by timefilter: %d of %d\n',...
    evbefore-evafter, evbefore)

if isempty(eventTimes)
    fprintf('eventTimes is empty\n');
    return
end

% Remove triggering events that are too close to the beginning or end
while eventTimes(1,1)<(epStartTime+win(1))
    eventTimes(1,:) = [];
end
while eventTimes(end,1)>(epEndTime-win(2))
    eventTimes(end,:) = [];
end

%% stack the event trig LFP
LFPtime=(epStartTime:1/samprate:epEndTime)';
evIdx = lookup(eventTimes(:,1),LFPtime);
for r=1:length(evIdx)
    for n = 1:length(nts)
        for l = 1:length(LFPtypes)
            lfpTypeIdx = find(contains(varargin(1:2:end), LFPtypes{l}));
            e = varargin{o(lfpTypeIdx)};
            out.data{l}{r}(n,:) = e{day}{ep}{nts(n)}.data(evIdx(r)-w(1):evIdx(r)+w(2));
        end
    end
end

%%
out.LFPtypes = LFPtypes;
out.LFPtimes = LFPtime;
out.eventStartIndices = evIdx;
out.epEndTime = epEndTime;
out.epStartTime = epStartTime;
out.clockrate = eeg{day}{ep}{idx(1,3)}.clockrate;
out.samprate = samprate;
% out.eventEndIndices = endidx;
out.win = win;
out.excludeperiods = timeFilter;
out.index = idx;
end

function out = init_out()
    out.data = [];
    out.LFPtypes = [];
    out.LFPtimes = [];
    out.eventStartIndices = [];
    out.epEndTime = [];
    out.epStartTime = [];
    out.clockrate = [];
    out.samprate = [];
%     out.eventEndIndices = [];
    out.win = [];
    out.excludeperiods = [];
    out.index = [];
end