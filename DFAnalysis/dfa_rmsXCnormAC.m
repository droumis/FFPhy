
function out = dfa_rmsXCnormAC(idx, excludeIntervals, varargin)
% get rms of xcnormac
% compatible with iterator: singledayanal 
% DR_2021-01-07

fprintf('day %d \n', idx(1));

requiredData = {'ca1rippleskons', 'lick', 'DIO', 'task', 'pos'};
varargins_allowed = ['animal', 'win', requiredData(:)'];

check_required(requiredData, varargin)
check_allowed(varargins_allowed, varargin)

if ~isempty(varargin)
    assign(varargin{:});
end

out = init_out(idx, animal);
day = idx(1);
eps = idx(2:end);

% get xp time series and rippwr timeseries
[xp, rip, pos, time] = collect_across_eps(lick, ca1rippleskons, pos, day, eps, ...
    excludeIntervals);
xpI = timestampsToIndicator(xp, time);

% get rewards well visits
% this should be replaced by well visit intervals
rewTimes = get_rewardTimes(task, DIO, day, eps);
intervalLen = 7; % s duration
rewIntervs = [rewTimes(:,1) rewTimes(:,1)+intervalLen];
if any((rewIntervs(2:end,1) > rewIntervs(1:end-1,2)))
    error(fprintf('overlapping intervals\n'))
end

[~, boutTimes] = getLickBoutLicks(animal, [repmat(day, length(eps),1) eps']);
boutTimesall = vertcat(boutTimes{day}{eps});

clc
plot(pos(:,1), pos(:,[2 3 4]))
for rI = 1:size(rewIntervs, 1)
    line([rewIntervs(rI,1), rewIntervs(rI,1)], ylim, 'color', 'g', 'linewidth', 4)
    line([rewIntervs(rI,2), rewIntervs(rI,2)], ylim, 'color', 'm', 'linewidth', 4)
end

% drop times not at reward
xpS = {};
ripS = {};
xc = {};
xcRMS = {};
for rI = 1:size(rewIntervs, 1)
    incT = isIncluded(time, rewIntervs(rI,:));
    xpS{1,rI} = xpI(incT);
    ripS{1,rI} = rip(incT);
    xc{rI} = xcNormAc(xpS{1,rI}, ripS{1,rI}, 'bin', 500, 'step', 50);
    xcRMS{rI} = rms(xc{rI});
end
% compute XCnormAC for all rewarded visits
d = nanmean(cell2mat(xc), 2);
m = nanmean(cell2mat(xcRMS));
%% compute RMS for 1000 permutations
% circularly permute within each well visit interval
perm_xc = {};
perm_xcRMS = [];

for i = 1:100
    for rI = 1:size(rewIntervs, 1)
        cutidx = randi([1,size(xpS{rI},1)], 1);
        pxpS = [xpS{1,rI}(cutidx:end); xpS{1,rI}(1:cutidx-1)];
        perm_xc{i,rI} = xcNormAc(pxpS, ripS{1,rI}, 'bin', 500, 'step', 50);
        perm_xcRMS(i, rI) = rms(perm_xc{i,rI});
    end
end
figure
histogram(perm_xcRMS(:),'numbins',100)

e = sum(perm_xcRMS < m);
line([m m], ylim)

% figure(1)
% plot(xcNormAc_full, 'color', 'k', 'linewidth', 4)
% hold on

out.xcRMS = xcRMS;
out.xc = xc;
out.perm_xc = perm_xc;
out.perm_xcRMS = perm_xcRMS;
end

function [xp, rip, posall, time] = collect_across_eps(lick, ca1rippleskons, ...
    pos, day, eps, excludeIntervals)
xp = [];
rip = [];
time = [];
posall = [];
for e = 1:eps
    ixp = lick{day}{eps(e)}.starttime;
    incxp = ~isIncluded(ixp, excludeIntervals);
    xp = [xp; ixp(incxp)];
    rip = [rip; ca1rippleskons{day}{eps(e)}{1}.powertrace'];
    time = [time; ca1rippleskons{day}{eps(e)}{1}.eegtimesvec_ref'];
    d = getfieldIdx(pos{day}{eps(e)}.fields, ...
        {'time' 'x-loess' 'y-loess' 'vel-loess'});
    posall = [posall; pos{day}{eps(e)}.data(:,d)];
end
end

function o = getfieldIdx(allfields, fields2get)
    f = split(allfields, ' '); 
    o = [];
for g = 1:length(fields2get)
    o = [o find(ismember(f, fields2get{g}))];
end
end

function out = init_out(idx, animal)
out.index = idx;
out.animal = animal;
out.xcNormAc = [];
end

