
function out = dfa_XCfromRew(idx, excludeIntervals, varargin)
% get xcnormac as a function of time since reward
% use with singledayanal
% DR_2020-12-30

fprintf('day %d \n', idx(1));

requiredData = {'ca1rippleskons', 'lick', 'DIO', 'task'};
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
[xp, rip, time] = collect_across_eps(lick, ca1rippleskons, day, eps, ...
    excludeIntervals);
% make binary time series of xp
xpidx = lookup(xp, time');
xpt = zeros(length(time),1);
xpt(xpidx)=1;

rewTimes = get_rewardTimes(task, DIO, day, eps);

intervalLen = 7; % s duration
rewIntervs = [rewTimes(:,1) rewTimes(:,1)+intervalLen];
if any((rewIntervs(2:end,1) - rewIntervs(1:end-1)) <= 0)
    error(fprintf('overlapping intervals\n'))
end
incT = (isIncluded(time, rewIntervs));
xpS = xpt(incT);   
ripS = rip(incT);
% cla
% compute XCnormAC for all rewarded visits
xcNormAc_full = xcNormAc(xpS, ripS, 'bin', 1500, 'step', 375);
% figure(1)
% plot(xcNormAc_full, 'color', 'k', 'linewidth', 4)
% hold on
%%
% partial exclusive xcNormAc
xcNormAc_incr = {};

for i = 1:intervalLen
    rewIntervs_incr = [rewTimes(:,1)+i-1 rewTimes(:,1)+i];
    incT = isIncluded(time, rewIntervs_incr);
    xpS = xpt(incT);   
    ripS = rip(incT);
    xcNormAc_incr{1,i} = xcNormAc(xpS, ripS, 'bin', 1500, 'step', 1500)';
%     plot(xcNormAc_incr{1,i}, 'color', rand(1,3), 'linewidth', 2)
end
% pause

out.xcNormAc_incr = xcNormAc_incr;
out.xcNormAc_full = xcNormAc_full;
end

function [xp, rip, time] = collect_across_eps(lick, ca1rippleskons, ...
    day, eps, excludeIntervals)
xp = [];
rip = [];
time = [];
for e = 1:eps
    ixp = lick{day}{eps(e)}.starttime;
    incxp = ~isIncluded(ixp, excludeIntervals);
    xp = [xp; ixp(incxp)];
    rip = [rip; ca1rippleskons{day}{eps(e)}{1}.powertrace'];
    time = [time; ca1rippleskons{day}{eps(e)}{1}.eegtimesvec_ref'];
end
end

function out = init_out(idx, animal)
out.index = idx;
out.animal = animal;
out.xcNormAc_incr = [];
out.xcNormAc_full = [];
end

