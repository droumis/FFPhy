function out = dfa_riptriglfp(idx, excludeperiods, varargin)

% gathers lfp around ripple times

win = [0.5 0.5];
appendindex = 1;
eventtype = 'rippleskons';
LFPtypes = [];
out.data = [];

if ~isempty(varargin)
    assign(varargin{:});
end

win = abs(win); % sometimes the preceding-event-window is negative..
if isempty(LFPtypes)
    disp('no lfp types given, using default: ''eeg'', ''ripple''')
    LFPtypes = {'eeg', 'ripple'};
end

% get the ripplekons events, which could be variably named
evid = find(contains(varargin(1:2:end), 'rippleskon'));
o = [1:2:length(varargin)]+1;
events = varargin{o(evid)};

eventtimes = [];
try
eventtimes = events{idx(1,1)}{idx(1,2)}{1}.starttime;
eventtimes = [eventtimes events{idx(1,1)}{idx(1,2)}{1}.endtime];
catch
    fprintf('no events detected for day%d ep%d\n', idx(1,1),idx(1,2))
    out = [];
    return
end

ecbefore = size(eventtimes,1);
eventtimes = eventtimes(~isExcluded(eventtimes(:,1),excludeperiods),:);
ecafter = size(eventtimes,1);
fprintf('%d of %d events discarded bc excluded periods in timefilter: d%d e%d\n',...
    ecbefore-ecafter, ecbefore, idx(1,1),idx(1,2))

if isempty(eventtimes)
    out.LFPtypes = LFPtypes;
    out.LFPtimes = [];
    out.eventStartIndices = [];
    out.eventEndIndices = [];
    out.win = win;
    out.index = idx;
    return
end

e = eeg{idx(1,1)}{idx(1,2)}{idx(1,3)}.data';
num_samp = length(e);
starttime = double(eeg{idx(1,1)}{idx(1,2)}{idx(1,3)}.starttime);
endtime_eeg = double(eeg{idx(1,1)}{idx(1,2)}{idx(1,3)}.endtime);
clockrate = eeg{idx(1,1)}{idx(1,2)}{idx(1,3)}.clockrate;
samprate = eeg{idx(1,1)}{idx(1,2)}{idx(1,3)}.samprate;
w = win*samprate;
endtime = (num_samp/samprate) + starttime;

%Remove triggering events that are too close to the beginning or end
while eventtimes(1,1)<(starttime+win(1))
    eventtimes(1,:) = [];
end
while eventtimes(end,2)>(endtime-win(2))
    eventtimes(end,:) = [];
end

LFPtimes=(starttime:1/samprate:endtime)';
stidx = lookup(eventtimes(:,1),LFPtimes);
endidx = lookup(eventtimes(:,2),LFPtimes);
%     windowStartIndices = lookup(eventtimes(:,1)-win(1),LFPtimes);
%     windowEndIndices = lookup(eventtimes(:,2)+win(2),LFPtimes);
%     eventStartLFPtimes = LFPtimes(eventStartIndices);
%     windowStartEndTimes = [eventStartLFPtimes(:)-win(1) eventStartLFPtimes(:)+win(2)];

%% STEP 3: Gather the ripple window data from all the regions
for r=1:length(stidx) %for each ripple from source region
    clear YripLFPdata
    for n = 1:length(idx(:,3))
        for l = 1:length(LFPtypes)
            out.data{l}{r}(n,:) = eval([LFPtypes{l}...
    '{idx(n,1)}{idx(n,2)}{idx(n,3)}.data(stidx(r)-w(1):stidx(r)+w(2));']);
        end
    end
end
out.LFPtypes = LFPtypes;
out.LFPtimes = LFPtimes;
out.eventStartIndices = stidx;
out.eventEndIndices = endidx;
out.win = win;
out.excludeperiods = excludeperiods;
out.index = idx;
end