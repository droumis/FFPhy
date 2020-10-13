function out = dfa_pctILB(idx, timeFilter, varargin)
% get pct swr iLB vs eLB
% use with singledayanal
% notebook: pctILB_2020_10_02.m
% DR_2020-10-02

fprintf('day %d \n',idx(1))
reqData = {'ca1rippleskons', 'lick'};
check_required(reqData, varargin)

stdthresh = [2 4 6 8];
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

maxthresh = maxthresh(incSWRs);
if isempty(swrTimes)
    return
end

%% Get licks in lick-burst intervals
[intraBoutXP, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
    varargin{:});
boutTimes = cell2mat({boutTimes{day}{eps}}');
% fprintf('%d lick bursts \n', size(boutTimes,1))

%% get iLB swrs

for istd = 1:length(stdthresh)
    stdidx = maxthresh > stdthresh(istd);
    swrTimes_std = swrTimes(stdidx);
    burstSWRTimes = swrTimes_std(isIncluded(swrTimes_std, boutTimes));
%     burstSWRTimes = sort(burstSWRTimes(~isnan(burstSWRTimes)));
    fprintf('::>%d Thresh:: %d of %d swrs within %d lickbouts. (%.02f pct swrs) \n', ...
        stdthresh(istd), numel(burstSWRTimes), numel(swrTimes_std), size(boutTimes,1),...
        length(burstSWRTimes)/length(swrTimes_std)*100)
%     if isempty(burstSWRTimes)
%         fprintf('no swrs in lick bouts for %d %d.. \n', day, ep)
%         keyboard
%     end
        
    out.numSWR(istd) = length(swrTimes_std);
    out.numSWRiLB(istd) = length(burstSWRTimes);
    out.burstSWRTimes{istd} = burstSWRTimes;
    out.pctSWRiLB(istd) = length(burstSWRTimes)/length(swrTimes_std)*100;
%     out.maxthresh = maxthresh;
end
end


function out = init_out()
out.index = [];
out.animal = [];
out.numSWR = [];
out.numSWRiLB = [];
out.burstSWRTimes = [];
out.pctSWRiLB = [];
% out.maxthresh = [];
end