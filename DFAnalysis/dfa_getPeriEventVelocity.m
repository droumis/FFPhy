%{
gathers velocity around event times for use with filterframework
Args:
    idx: (array) day epoch index [day, epoch]
    excludeperiods: (array) periods in epoch to exclude [start1, end1; ...]
    events: (cellarray) {day}{epoch}{1}.starttime .endtime
    position: (cellarray) {day}{epoch}.data .fields

compatible with iterator: singleepochanal.m

Author: Demetris Roumis 2019
%}

function out = dfa_getPeriEventVelocity(idx, excludeperiods, events, pos, varargin)

win = [0.5 0.5];
appendindex = 1;
eventtype = 'rippleskons';
if ~isempty(varargin)
    assign(varargin{:});
end

out.data = [];
win = abs(win); % assert positive window duration
day = idx(1,1);
epoch = idx(1,2);

out.velTimes = [];
out.eventStartIdx = [];
out.eventEndIdx = [];
out.win = win;
out.index = idx;

% check for events
eventTime = [];
try
    eventTime = [events{day}{epoch}{1}.starttime events{day}{epoch}{1}.endtime];
catch
    fprintf('no events detected for day%d ep%d\n', day,epoch)
    return
end
% print proportion of included events 
ecbefore = size(eventTime,1);
eventTime = eventTime(~isExcluded(eventTime(:,1),excludeperiods),:);
ecafter = size(eventTime,1);
fprintf('%d / %d events included : d%d e%d\n', ecafter, ecbefore, day,epoch)
if isempty(eventTime)
    return
end

% p = pos{day}{epoch}.data;
timecol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'time'));
% velcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'vel-loess'));
% xcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'x-loess'));
% ycol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'y-loess'));
% dircol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'dir-loess'));
% timevel = pos{day}{epoch}.data(:, [timecol velcol xcol ycol dircol]);

% num_samp = size(pos{day}{epoch}.data,1);
starttime = pos{day}{epoch}.data(1,timecol);
endtime = pos{day}{epoch}.data(end,timecol);
% clockrate = eeg{day}{epoch}{idx(1,3)}.clockrate;
samprate = round(1/(mode(diff(pos{day}{epoch}.data(:,timecol))))); % Hz
w = win*samprate;
% endtime = (num_samp/samprate) + starttime;

%Remove events that are too close to the beginning or end
while eventTime(1,1)<(starttime+win(1))
    eventTime(1,:) = [];
end
while eventTime(end,2)>(endtime-win(2))
    eventTime(end,:) = [];
end

% epochtime=pos{day}{epoch}.data(:,timecol);
% winStartIdx = lookup(eventtimes(:,1)-win(1),pos{day}{epoch}.data(:,timecol));
% winEndIdx = lookup(eventtimes(:,2)+win(2),pos{day}{epoch}.data(:,timecol));
stidx = lookup(eventTime(:,1),pos{day}{epoch}.data(:,timecol));
endidx = lookup(eventTime(:,2),pos{day}{epoch}.data(:,timecol));
%     eventStartLFPtimes = LFPtimes(eventStartIndices);
%     windowStartEndTimes = [eventStartLFPtimes(:)-win(1) eventStartLFPtimes(:)+win(2)];
out.data = zeros(size(pos{day}{epoch}.data,2),sum(w)+1, length(stidx));
%% Gather the event window data
for r=1:length(stidx)
    out.data(:,:,r) = pos{day}{epoch}.data(stidx(r)-w(1):stidx(r)+w(2), :)';
end
out.field = strsplit(pos{day}{epoch}.fields,' ');
out.dims = {'field', 'time', 'event'};
out.time = (-w(1):w(2))*(1/samprate);
out.eventStartTimes = stidx;
out.eventEndTimes = endidx;
out.win = win;
out.excludeperiods = excludeperiods;
out.index = idx;
end