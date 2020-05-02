function out = dfa_phaseXcorr(idx, timeFilter, varargin)
% filterframework analysis function 
% return pairwise phase cross-correlation
% Options:
%   'bin', n 	binsize in sec. Default 0.002 (2 ms)
%   'tmax', n	maximum time for cross correlation. Default 1 sec.

bin = .001; % s
tmax = .5; % s
try assign(varargin{:}); catch; end

% check for required data in varargin
reqData = {'spikes', 'licks'};
for s = 1:length(reqData)
    if ~any(cell2mat(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 0)))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end

% fuck this won't work for pairs.. 
day = idx(1);
eps = idx(4:end);
nt = idx(2);
clust = idx(3);

% for each cell we calculate the cross correlation 
try
    t1 = spikes{idx(1)}{idx(2)}{idx(3)}{idx(4)}.data;
    t2 = spikes{idx(1)}{idx(2)}{idx(5)}{idx(6)}.data;
catch
    % if either of those produced an error, we return NaNs 
    return;
end
%apply the exclude rule
t1inc = [];
t2inc = [];
if (length(t1))
    t1inc = t1(find(~isExcluded(t1(:,1), timeFilter)),1);
else
    return
end
if (length(t2))
    t2inc = t2(find(~isExcluded(t2(:,1), timeFilter)),1);
else
    return
end

% %% get spikes
% spikeTimes = [];
% numSpikesPerEp = [];
% for e = 1:length(eps)
%     try
%         spikeTimes = [spikeTimes; spikes{day}{eps(e)}{nt}{clust}.data(:,1)];
%         numSpikesPerEp = [numSpikesPerEp size(spikes{day}{eps(e)}{nt}{clust}.data,1)];
%     catch
%         continue
%     end
% end



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


out.xcc1vc2

end

function out = init_out()
out.tmax = [];
out.bin = [];

end