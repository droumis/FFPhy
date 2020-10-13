
function out = dfa_XPtrigAvgRip(idx, timeFilter, varargin)
% get XP trig avg of rip power
% use with singledayanal
% DR_2020-09-30

fprintf('day %d \n',idx(1))
reqData = {'ca1rippleskons', 'lick'};
check_required(reqData, varargin)

win = [-.5 .5];
eventType = 'ca1rippleskons';
srate = 1500;
if ~isempty(varargin)
    assign(varargin{:});
end

out = init_out(); % init output
out.index = idx;
out.animal = animal;
day = idx(1);
eps = idx(2:end);

%% get rip pwr trace
evid = find(contains(varargin(1:2:end), eventType));
o = [1:2:length(varargin)]+1;
swr = varargin{o(evid)};

rippwr = [];
time = [];

for e = 1:eps
    rippwr = [rippwr; swr{day}{eps(e)}{1}.powertrace'];
    time = [time; swr{day}{eps(e)}{1}.eegtimesvec_ref'];
end

% zscore
rippwr = rippwr';
zrippwr = ((rippwr - nanmean(rippwr))./nanstd(rippwr));

%% Get licks in lick-burst intervals
[intraBoutXP, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
    varargin{:});
boutTimes = cell2mat({boutTimes{day}{eps}}');
intraBoutXP = cell2mat({intraBoutXP{day}{eps}}');
intraBoutXP = intraBoutXP(:,1);
fprintf('%d XP within %d bursts \n', numel(intraBoutXP), size(boutTimes,1))

% for each iLB XP, get a zrippwr trace within window
% return this stack
% also return the nanmean and nanstd of the stack
zrippwr_XPtrig = [];
rippwr_XPtrig = [];
for iXP = 1:length(intraBoutXP)
    xptimeIdx = min(find(time>=intraBoutXP(iXP)));
    if diff([intraBoutXP(iXP) time(xptimeIdx)]) > .001
        sprintf('%.03f diff, skipping',diff([intraBoutXP(iXP) time(xptimeIdx)]))
        continue
    end
    startIdx = xptimeIdx-abs(win(1)*srate)-1;
    endIdx = xptimeIdx+abs(win(2)*srate);
    try
        zrippwr_XPtrig = [zrippwr_XPtrig; zrippwr(startIdx:endIdx)];
        rippwr_XPtrig = [rippwr_XPtrig; rippwr(startIdx:endIdx)];
    catch
        disp('window exceeds bounds, skipping')
        continue
    end
end

out.zrippwr_XPtrig = zrippwr_XPtrig;
out.mean_zrippwr_XPtrig = nanmean(zrippwr_XPtrig,1);
out.sem_zrippwr_XPtrig = nanstd(zrippwr_XPtrig,1)/sqrt(size(zrippwr_XPtrig,1));
out.rippwr_XPtrig = rippwr_XPtrig;
out.mean_rippwr_XPtrig = nanmean(rippwr_XPtrig,1);
out.sem_rippwr_XPtrig = nanstd(rippwr_XPtrig,1)/sqrt(size(rippwr_XPtrig,1));
out.win = win;
out.srate = srate;
out.time = linspace(win(1),win(2),length(out.mean_zrippwr_XPtrig));
end

function out = init_out()
out.index = [];
out.animal = [];
out.time = [];

out.zrippwr_XPtrig = [];
out.mean_zrippwr_XPtrig = [];
out.sem_zrippwr_XPtrig = [];
out.rippwr_XPtrig = [];
out.mean_rippwr_XPtrig = [];
out.sem_rippwr_XPtrig = [];
end