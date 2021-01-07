

function out = dfa_XPtrigSWR(idx, timeFilter, varargin)


% get rip position
% use with singledayanal
% DR_2020-10-12

fprintf('day %d \n',idx(1))
reqData = {'ca1rippleskons'};
check_required(reqData, varargin)

win = [-1 1];
bin = .001;
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
SWRTimes = [];
XP = [];
maxthresh = [];
swrEnd = [];
swr_ep = [];
for e = 1:length(eps)
    ep = eps(e);
    try
        SWRTimes = [SWRTimes; swr{day}{ep}{1}.starttime];
        XP = [XP; lick{day}{ep}.starttime];
        swr_ep = [swr_ep; repmat(ep,length(SWRTimes),1)];
        maxthresh = [maxthresh; swr{day}{ep}{1}.maxthresh];
        swrEnd = [swrEnd; swr{day}{ep}{1}.endtime];
    catch
        try
            SWRTimes = [SWRTimes; swr{day}{ep}.starttime];
            XP = [XP; lick{day}{ep}.starttime];
            swr_ep = [swr_ep; repmat(ep,length(SWRTimes),1)];
            maxthresh = [maxthresh; swr{day}{ep}.maxthresh];
            swrEnd = [swrEnd; swr{day}{ep}.endtime];
        catch
            fprintf('no swr events detected for day%d ep%d\n', day, ep)
            return
        end
    end
end
% apply timefilter to swrs
incSWRs = ~isIncluded(SWRTimes(:,1), timeFilter);
fprintf('%d of %d swr events discarded bc excluded periods in timefilter: d%d\n',...
    sum(incSWRs), length(incSWRs), day)
SWRTimes = SWRTimes(incSWRs,:);
SWRTimes_maxThresh = maxthresh(incSWRs,:);
swr_ep = swr_ep(incSWRs);

if isempty(SWRTimes)
    return
end

%% Get XP
[XP_iLB, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
    varargin{:});
XP_iLB = cell2mat({XP_iLB{day}{eps}}');
XP_iLB = XP_iLB(:,1);
%% Get SWR offset from nearest XP
XPnearSWR_idx = knnsearch(XP_iLB, SWRTimes);
SWRtimeFromNearestXP = SWRTimes - XP_iLB(XPnearSWR_idx);

r = rand(length(SWRTimes),1);
rSt = -.5;
rEnd = .5;
rjitt = (rEnd-rSt) * r + rSt;
SWRTimes_Sh = SWRTimes + rjitt;
XPnearSWR_idx_Sh = knnsearch(XP_iLB, SWRTimes_Sh);
SWRtimeFromNearestXP_Sh = SWRTimes_Sh - XP_iLB(XPnearSWR_idx_Sh);
%% get window of SWR events around each XP
time = -abs(win(1))-(bin/2) : bin : win(2)+(bin/2);
XPiLBtrigSWR_raster = cell2mat(arrayfun(@(r) histc(SWRTimes, r + time), ...
    XP_iLB, 'un', 0)')';

XPtrigSWR_raster = cell2mat(arrayfun(@(r) histc(SWRTimes, r + time), ...
    XP, 'un', 0)')';
%% output
out.time = time;
out.SWRtimeFromNearestXP = SWRtimeFromNearestXP;
out.SWRtimeFromNearestXP_Sh = SWRtimeFromNearestXP_Sh;
out.XPiLBtrigSWR_raster = XPiLBtrigSWR_raster;
out.XPtrigSWR_raster = XPtrigSWR_raster;
end

function out = init_out()
out.index = [];
out.animal = [];

out.time = [];
out.SWRtimeFromNearestXP = [];
out.SWRtimeFromNearestXP_Sh = [];
out.XPiLBtrigSWR_raster = [];
out.XPtrigSWR_raster = [];
end