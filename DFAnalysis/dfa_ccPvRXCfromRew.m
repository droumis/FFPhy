function out = dfa_ccPvRXCfromRew(idx, excludeIntervals, varargin)
% calculate correlation coefficient between XCnormAc-predicted ripple 
% power trace vs the real, as a function of time since reward 
%
% use with singledayanal
%
%
% DKR 2021-01-01

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

incT = isIncluded(time, rewIntervs);
xpS = xpt(incT);   
ripS = rip(incT);
% compute XCnormAC for all rewarded visits
xcnorm = xcNormAc(xpS, ripS, 'bin', 1500, 'step', 250);

% Convolve the xcnorm with the BB pulse signal
ripPredicted = conv(xpt, xcnorm, 'same');

% get the corrcoef for each time bin of each rewarded visit
coefs = nan(size(rewIntervs,1), intervalLen);
for ibin = 1:intervalLen
    for ivisit = 1:size(rewIntervs,1)
        istart = rewIntervs(ivisit, 1)+ibin-1;
        iend = rewIntervs(ivisit, 1)+ibin;
        incT = isIncluded(time, [istart iend]);
        predS = ripPredicted(incT);
        ripS = rip(incT);
        predS = predS - nanmean(predS);
        ripS = ripS - nanmean(ripS);
        [coef, pval] = corr(predS, ripS);
        coefs(ivisit, ibin) = coef;
    end
end
out.coefs = coefs;

rewIdx = lookup(rewIntervs(:,1), time);
t = 0:1500*intervalLen-1;

out.predRipStack = ripPredicted(rewIdx+t);
out.realRipStack = rip(rewIdx+t);
out.rewIntervs = rewIntervs;
out.xcnorm = xcnorm;
out.rip = rip;
out.xp = xp;
out.time = time;

end


function [xp, rip, time] = collect_across_eps(lick, ca1rippleskons, ...
    day, eps, excludeIntervals)
xp = [];
rip = [];
time = [];
for e = 1:length(eps)
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


end