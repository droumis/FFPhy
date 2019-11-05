function out = dfa_eventTrigLFP(idx, excludeIntervals, varargin)

% gathers lfp around event times
% DR 19

fprintf('%d %d %d %d\n',idx)
win = [0.5 0.5];
eventType = 'lick';
LFPtypes = [];
if ~isempty(varargin)
    assign(varargin{:});
end

day = idx(1,1);
ep = idx(1,2);
nts = idx(:,3);
out = init_out();
win = abs(win); % sometimes the preceding-event-window is negative..

%% LFP

if isempty(LFPtypes)
    disp('no lfp types given, using default: ''eeg'', ''ripple''')
    LFPtypes = {'eeg', 'ripple'};
end
lfpTypeIdx = find(contains(varargin(1:2:end), LFPtypes{1}));
o = [1:2:length(varargin)]+1;
e = varargin{o(lfpTypeIdx)};
num_samp = length(e{day}{ep}{nts(end)}.data);
epStartTime = double(e{day}{ep}{nts(end)}.starttime);
samprate = e{day}{ep}{nts(end)}.samprate;
epEndTime = (num_samp/samprate) + epStartTime;
w = win*samprate;

%% get events from ~excludeperiods
evid = find(contains(varargin(1:2:end), eventType));
events = varargin{o(evid)};
eventtimes = [];
try
    eventtimes = events{day}{ep}{1}.starttime;
catch
    try
        eventtimes = events{day}{ep}.starttime;
    catch
        fprintf('no events detected for day%d ep%d\n', day,ep)
        return
    end
end

ecbefore = size(eventtimes,1);
eventtimes = eventtimes(~isExcluded(eventtimes(:,1),excludeIntervals),:);
ecafter = size(eventtimes,1);
fprintf('%d of %d events discarded bc excluded periods in timefilter: d%d e%d\n',...
    ecbefore-ecafter, ecbefore, day,ep)

if isempty(eventtimes)
    return
end

% Remove triggering events that are too close to the beginning or end
while eventtimes(1,1)<(epStartTime+win(1))
    eventtimes(1,:) = [];
end
while eventtimes(end,1)>(epEndTime-win(2))
    eventtimes(end,:) = [];
end

%% stack the event trig LFP
LFPtime=(epStartTime:1/samprate:epEndTime)';
evIdx = lookup(eventtimes(:,1),LFPtime);
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
out.excludeperiods = excludeIntervals;
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