

function out = dfa_XPprepostSWR(idx, timeFilter, varargin)
% get XP offset, pre and post SWR
% use with singledayanal
% DR_2020-09-30

fprintf('day %d \n',idx(1))
reqData = {'ca1rippleskons', 'lick'};
check_required(reqData, varargin)

run_shuffle = 1;
eventType = 'ca1rippleskons';
if ~isempty(varargin)
    assign(varargin{:});
end

out = init_out(); % init output
out.index = idx;
out.animal = animal;
day = idx(1);
eps = idx(2:end);
%% Get SWRs in timeFilter

evid = find(contains(varargin(1:2:end), eventType));
o = [1:2:length(varargin)]+1;
swr = varargin{o(evid)};
swrTimes = [];
maxthresh = [];
swrEnd = [];
for e = 1:length(eps)
    ep = eps(e);
    try
        swrTimes = [swrTimes; swr{day}{ep}{1}.starttime];
        maxthresh = [maxthresh; swr{day}{ep}{1}.maxthresh];
        swrEnd = [swrEnd; swr{day}{ep}{1}.endtime];
    catch
        try
            swrTimes = [swrTimes; swr{day}{ep}.starttime];
            maxthresh = [maxthresh; swr{day}{ep}.maxthresh];
            swrEnd = [swrEnd; swr{day}{ep}.endtime];
        catch
            fprintf('no swr events detected for day%d ep%d\n', day, ep)
            return
        end
    end
end
% apply timefilter to swrs
incSWRs = ~isIncluded(swrTimes(:,1), timeFilter);
fprintf('%d of %d swr events discarded bc excluded periods in timefilter: d%d\n',...
    sum(incSWRs), length(incSWRs), day)
swrTimes = swrTimes(incSWRs,:);
% out.maxthresh = maxthresh;
% out.swrEnd = swrEnd;

if isempty(swrTimes)
    return
end

%% Get licks in lick-burst intervals
[intraBoutXP, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
    varargin{:});
boutTimes = cell2mat({boutTimes{day}{eps}}');
intraBoutXP = cell2mat({intraBoutXP{day}{eps}}');
intraBoutXP = intraBoutXP(:,1);
fprintf('%d XP within %d bursts \n', numel(intraBoutXP), size(boutTimes,1))


%% get xp pre post swr
% get swr events within lick-burst intervals
burstSWRTimes = swrTimes(isIncluded(swrTimes, boutTimes));
burstSWRTimes = sort(burstSWRTimes(~isnan(burstSWRTimes)));
fprintf(' %d swrs within %d lickbouts. (%.02f pct swrs) \n', size(boutTimes,1), ...
    numel(burstSWRTimes), length(burstSWRTimes)/length(swrTimes)*100)
if isempty(burstSWRTimes)
    fprintf('no swrs in lick bouts for %d %d.. skipping\n', day, ep)
    return
end
% Get the swr-containing licks
% histc given licks as edges: get bin idx. lickedges(binidx) is the prior edge of bin, 
% aka the swr-preceding lick 
% [~,~,swrLickBinidx] = histcounts(burstSWRTimes, intraBoutXP);
% 
% % time since last lick
% swrTimeSinceLick = burstSWRTimes - intraBoutXP(swrLickBinidx);
% if any(swrTimeSinceLick > .250)
%     error('wut')
% end
% swr-containing ILI
for b = 1:length(burstSWRTimes)
    XPpreSWR(b,1) = max(intraBoutXP(intraBoutXP<burstSWRTimes(b)));
    XPpostSWR(b,1) = min(intraBoutXP(intraBoutXP>=burstSWRTimes(b)));
end
XPpreSWR_offset = burstSWRTimes - XPpreSWR;
XPpostSWR_offset = burstSWRTimes - XPpostSWR;

if run_shuffle
    r = rand(length(XPpreSWR),1);
    ILI = XPpostSWR - XPpreSWR; % ILI
    SWRtimes_shuf = XPpreSWR + (ILI .* r);% rand fraction of each ILI as shuf swr
    XPpreSWR_offset_shuf = SWRtimes_shuf - XPpreSWR;
    XPpostSWR_offset_shuf = SWRtimes_shuf - XPpostSWR;
end

out.XPpreSWR = XPpreSWR;
out.XPpostSWR = XPpostSWR;
out.XPpreSWR_offset = XPpreSWR_offset;
out.XPpostSWR_offset = XPpostSWR_offset;
out.XPpreSWR_offset_shuf = XPpreSWR_offset_shuf;
out.XPpostSWR_offset_shuf = XPpostSWR_offset_shuf;
out.SWRTimes_iLB = burstSWRTimes;
out.XPTimes_iLB = intraBoutXP;
end

function out = init_out()
out.index = [];
out.animal = [];

out.XPpreSWR = [];
out.XPpostSWR = [];
out.XPpreSWR_offset = [];
out.XPpostSWR_offset = [];
out.XPpreSWR_offset_shuf = [];
out.XPpostSWR_offset_shuf = [];
out.SWRTimes_iLB = [];
out.XPTimes_iLB = [];
end