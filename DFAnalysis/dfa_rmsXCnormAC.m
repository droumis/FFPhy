
function out = dfa_rmsXCnormAC(idx, excludeIntervals, varargin)
% get rms of xcnormac
% compatible with iterator: singledayanal 
% DR_2021-01-07

fprintf('day %d \n', idx(1));
bin = 500;
step = 150;

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

% get ing Intervals;
ingIntervals = getIngestionIntervals(pos, DIO, task, day, eps);
ingIntervals = vertcat(ingIntervals{:});

% get xp time series and rippwr timeseries
[xp, rip, pos, time] = collect_across_eps(lick, ca1rippleskons, pos, day, eps, ...
    excludeIntervals);
xpI = timestampsToIndicator(xp, time);

%%
xpS = {};
ripS = {};
xc = {};
xcRMS = {};
mid = round(bin / 2);
hwin = 100;
win4RMS = mid-hwin:mid+hwin-1;
for rI = 1:size(ingIntervals, 1)
    incT = isIncluded(time, ingIntervals(rI,:));
    xpS{1,rI} = xpI(incT);
    ripS{1,rI} = rip(incT);
    figure(1)
    ax(1) = subaxis(2,1,1);
    plot(time(incT), xpS{1,rI})
    ax(2) = subaxis(2,1,2);
    plot(time(incT), ripS{1,rI})
    linkaxes(ax, 'x')
    axis tight
    xc{rI} = xcNormAc(xpS{1,rI}, ripS{1,rI}, 'bin', bin, 'step', step);
    xcRMS{rI} = rms(xc{rI}(win4RMS));
end
% compute XCnormAC for all rewarded visits
d = nanmean(cell2mat(xc), 2);
m = nanmean(cell2mat(xcRMS));

% figure
% plot(d)
% compute RMS for 1000 permutations
% circularly permute within each well visit interval
perm_xc = {};
perm_xcRMS = [];
tic
for i = 1:100
    for rI = 1:size(ingIntervals, 1)
        cutidx = randi([1,size(xpS{rI},1)], 1);
        pxpS = [xpS{1,rI}(cutidx:end); xpS{1,rI}(1:cutidx-1)];
        perm_xc{i,rI} = xcNormAc(pxpS, ripS{1,rI}, 'bin', bin, 'step', step);
        perm_xcRMS(i, rI) = rms(perm_xc{i,rI}(win4RMS));
    end
end
toc
figure
histogram(perm_xcRMS(:),'numbins',100)
line([m m], ylim)
%%
e = sum(sum(perm_xcRMS < m));


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
    d = getFieldIndex(pos{day}{eps(e)}.fields, ...
        {'time' 'x-loess' 'y-loess' 'vel-loess'});
    posall = [posall; pos{day}{eps(e)}.data(:,d)];
end
end

function out = init_out(idx, animal)
out.index = idx;
out.animal = animal;
out.xcNormAc = [];
end

