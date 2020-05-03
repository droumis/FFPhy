function out = dfa_phaseXcorr(idx, timeFilter, varargin)
% filterframework analysis function 
% return pairwise phase cross-correlation
%
% Input:
%   idx : [day ntA clustA ntB clustB ep1..epN]
%
% Options:
%   'bin', n 	binsize in sec. Default 0.002 (2 ms)
%   'tmax', n	maximum time for cross correlation. Default 1 sec.

bin = .001; % s
tmax = .5; % s
try assign(varargin{:}); catch; end

% check for required data in varargin
reqData = {'spikes', 'lick'};
for s = 1:length(reqData)
    if ~any(cell2mat(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 0)))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end

% fuck this won't work for pairs.. 
day = idx(1);
ntA = idx(2);
clustA = idx(3);
ntB = idx(4);
clustB = idx(5);
eps = idx(6:end);

out = init_out();
% for each cell we calculate the cross correlation 
t1 = [];
t2 = [];
for iep = 1:length(eps)
    try
        t1 = [t1; spikes{day}{eps(iep)}{ntA}{clustA}.data];
        t2 = [t2; spikes{day}{eps(iep)}{ntB}{clustB}.data];
    catch
        warning('couldnt load spikes %s d%d [ntA clA ntB clB]:: %d %d %d %d\n', ...
            animal, day, nyA, clustA, ntB, clustB);
        return
        % if either of those produced an error, we return NaNs 
    end
end
%apply the exclude rule
t1inc = [];
t2inc = [];
timeFilterCat = cell2mat([{timeFilter{:}}']); % combine epoch timefilters
t1cnd = t1(~logical(isExcluded(t1, timeFilterCat)));
t2cnd = t2(~logical(isExcluded(t2, timeFilterCat)));

%% get lickbout licks
[intraBoutXP, ~] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], varargin);
intraBoutXPvec = intraBoutXP{day}{eps};
eventTimes = intraBoutXPvec(ismember(intraBoutXPvec, eventTimes));

numEventsPerEp = [];
for e = 1:length(eps)
    epStartTime = spikes{day}{eps(e)}{nt}{clust}.timerange(1);
    epEndTime = spikes{day}{eps(e)}{nt}{clust}.timerange(2);
    numEventsPerEp = [numEventsPerEp; sum(logical(isExcluded(eventTimes, [epStartTime epEndTime])))];
end

%% transform spikes into a phase vector
% how? each spike is a spike time.. i need each spike time relative to
% preceeding lick then ratio to interlickinterval duration to get %.. then
% radians.. then that + (2pi*[preceding lick number])



%% run phase xcorr

% get xc
xc = spikexcorr(t1cnd,t2cnd, bin, tmax);
% normalize
normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2);
% smooth
try
    nstd=round(sw1/(xc.time(2) - xc.time(1)));
catch
    if xc.nspikes1 == 0
        emptynt = index(3:4);
    elseif xc.nspikes2 == 0
        emptynt = index(5:6);
    end
    fprintf('no spikes %d %d %d %d\n', day, epoch, emptynt(1),emptynt(2))
    out = make_blank(index);
    return
end
g1 = gaussian(nstd, 2*nstd+1);
smthxc = smoothvect(normxc, g1);
% Get Nevents in raw xc from -200ms to 200ms
% bins = find(abs(xc.time)<=0.2);
%     Neventscorr = sum(xc.c1vsc2(bins));
% Expected probability
p1 = xc.nspikes1/T;
p2 = xc.nspikes2/T; % fr in Hz
expProb = p1*p2; % per sec


% compute the excess correlation at 0 lag
exc = excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, sw1, sw2);
% compute RMS
xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts);
% output
out.index = index;
% out.T = T;
out.xc = xc;
out.normxc = normxc;
out.smthxc = smthxc;
out.expProb = expProb;
% out.xcShfmean = xcShfmean;
% out.xcShfLag0 = xcShfLag0;
out.excesscorr = nanmean(exc); % nanmean in case there are two bins equally near lag zero
out.xcrms = xcrms;
end

function out = init_out()
out.tmax = [];
out.bin = [];
out.xc = [];
out.normxc = [];
out.smthxc = [];
out.expProb = [];
out.excesscorr = [];
out.xcrms = [];
end