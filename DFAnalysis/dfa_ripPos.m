function out = dfa_ripPos(idx, timeFilter, varargin)

% get rip position
% use with singledayanal
% DR_2020-10-12

fprintf('day %d \n',idx(1))
reqData = {'ca1rippleskons', 'pos','lick'};
check_required(reqData, varargin)

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
maxthresh = [];
swrEnd = [];
swr_ep = [];
for e = 1:length(eps)
    ep = eps(e);
    try
        SWRTimes = [SWRTimes; swr{day}{ep}{1}.starttime];
        swr_ep = [swr_ep; repmat(ep,length(SWRTimes),1)];
        maxthresh = [maxthresh; swr{day}{ep}{1}.maxthresh];
        swrEnd = [swrEnd; swr{day}{ep}{1}.endtime];
    catch
        try
            SWRTimes = [SWRTimes; swr{day}{ep}.starttime];
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

%% Get pos
posdata = [];
for e = 1:length(eps)
    ep = eps(e);
    posdata = [posdata; pos{day}{eps(e)}.data];
    posfields = pos{day}{eps(e)}.fields;
end
xstring = 'x-loess';
xcol = find(cell2mat(cellfun(@(x) strcmp(x,xstring), strsplit(posfields, ' '), ...
    'UniformOutput', false)));
ystring = 'y-loess';
ycol = find(cell2mat(cellfun(@(x) strcmp(x,ystring), strsplit(posfields, ' '), ...
    'UniformOutput', false)));

swr_posidx = knnsearch(posdata(:,1), SWRTimes(:,1));
SWRTimes_xpos = posdata(swr_posidx,xcol);
SWRTimes_ypos = posdata(swr_posidx,ycol);

%% Get lick-burst intervals
[intraBoutXP, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
    varargin{:});
boutTimes = cell2mat({boutTimes{day}{eps}}');

%% get swr events within lick-burst intervals
SWRTimes_iLB = isIncluded(SWRTimes, boutTimes);

%% output
out.txy = posdata(:,[1 xcol ycol]);
out.SWR_time = SWRTimes;
out.SWR_iLB = SWRTimes_iLB;
out.SWR_xpos = SWRTimes_xpos;
out.SWR_ypos = SWRTimes_ypos;
out.SWR_maxThresh = SWRTimes_maxThresh;
end

function out = init_out()
out.index = [];
out.animal = [];

out.txy = [];
out.SWR_time = [];
out.SWR_iLB = [];
out.SWR_xpos = [];
out.SWR_ypos = [];
out.SWR_maxThresh = [];
end