function out = dfa_rewTrigSWRXP(idx, timeFilter, varargin)

% - reward triggered swr and licks
% - for use with singledayanal
% DR_2020-10-13

fprintf('day %d \n',idx(1))
reqData = {'ca1rippleskons', 'DIO','lick', 'task'};
check_required(reqData, varargin)
win = [-2 10];
bin = .001;
eventType = 'ca1rippleskons';
if ~isempty(varargin)
    assign(varargin{:});
end

out = init_out();
out.index = idx;
out.animal = animal;
day = idx(1);
eps = idx(2:end);

%% get DIO reward output

[dioOut, dioOutfields] = get_dio_output(DIO, task, day, eps);

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
% [intraBoutXP, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
%     varargin{:});
% boutTimes = cell2mat({boutTimes{day}{eps}}');
% intraBoutXP = cell2mat({intraBoutXP{day}{eps}}');
% intraBoutXP = intraBoutXP(:,1);
% fprintf('%d XP within %d bursts \n', numel(intraBoutXP), size(boutTimes,1))
xp = [];
for e = 1:length(eps)
    ep = eps(e);
    try
        xp = [xp; lick{day}{ep}.eventtime];
    catch
        xp = [xp; lick{day}{ep}.starttime]; % legacy name
    end
end

time = -abs(win(1))-(bin/2) : bin : win(2)+(bin/2);

swr_prth = cell2mat(arrayfun(@(r) histc(swrTimes, r + time), dioOut(:,5), 'un', 0)')';
xp_prth = cell2mat(arrayfun(@(r) histc(xp, r + time), dioOut(:,5), 'un', 0)')';

%% output
out.rewTrigSWR = swr_prth;
out.rewTrigXP = xp_prth;
out.time = time;

end

function out = init_out()
out.index = [];
out.animal = [];

out.rewTrigSWR = [];
out.rewTrigXP = [];
out.time = [];
end

