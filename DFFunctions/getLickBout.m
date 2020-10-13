

function out = getLickBout(datadir, animal, epochs, varargin)
% out = getLickBout(datadir, animal, epochs, varargin)
% Get the lick bout intervals, and time from last lick.
% This is a timefilter function, designed for FilterFramework.
% inputs:
% - datadir
% - animal
% - epochs = [day epoch; ...]
%
% output is a day, epoch-nested struct of time vec and var vecs
% var vecs are:
% - timeFromLick (continuous, seconds)
% - lickBout (binary, indicating inclusion in a lick bout)

%{
Notes:
get_licks creates the lick data structure that this uses

'datadir' arg1 currently needs to be there bc of the way setfiltertime
calls the time funcs... but i don't use it so outside of setfiltertime i
should just pass an empty arg there..

@DKR
%}
output_intervals = 0;
maxILIthresh = 1; % max burst ili threshold in seconds
minILIthresh = .06; % min burst ili threshold in seconds
minBoutLicks = 2; %filter out bouts with less than boutNum licks
lick = [];
if ~isempty(varargin)
    assign(varargin{:})
end

if isempty(lick)
    andef = animaldef(animal);
    days = unique(epochs(:,1));
    lick = loaddatastruct(andef{2}, animal, 'lick', days);
end

out = {};
for i = 1:size(epochs,1)
    day = epochs(i,1);
    epoch = epochs(i,2);
    if isa(lick, 'cell')
        % get lick bouts
        try
            lickTime = lick{day}{epoch}.eventtime;
        catch
            lickTime = lick{day}{epoch}.starttime; % legacy
        end
    else
        lickTime = lick;
    end
    % keep valid lickburst licks
    while 1
        lickTimesILI = diff(lickTime);
        if any(lickTimesILI < minILIthresh)
            lickTime(find(lickTimesILI < minILIthresh)+1) = [];
        else
           break
        end
    end
    % get burst intervals
    ili = diff(lickTime);
    g = [ili < maxILIthresh];
    bi = vec2list(g, 1:length(g));
    boutIntvIncl = (diff(bi,[],2)+1) >= minBoutLicks;
    bii = bi(boutIntvIncl,:);
    % convert ili indices into time
%     o = [lickTime(bii(:,1)), lickTime(bii(:,2)+1)];
    boutIntvStart = lickTime(bii(:,1));
    boutIntvEnd = lickTime(bii(:,2)+1);
    lickBoutXP = lickTime(isIncluded(lickTime, [boutIntvStart boutIntvEnd]));
%     % 
%     boutIntvStart = lickTime(find(diff([ili < maxILIthresh]) == 1)+1);
%     boutIntvEnd = lickTime(find(diff([ili > maxILIthresh]) == 1)+1);
%     while boutIntvStart(1) > boutIntvEnd(1)
%         boutIntvEnd(1) = [];
%         if boutIntvEnd(end)<boutIntvStart(end)
%             boutIntvStart(end) = [];
%         end
%     end
%     % filter out bouts with less than minBoutLicks
%     assert(length(boutIntvStart)==length(boutIntvEnd),'bout start end dims are unequal')
%     try
%         licksInbout = isIncluded(lickTime, boutIntv);
%     catch
%         fprintf('error defining lick bouts for %d %d\n', day, epoch)
%         continue
%     end
%     N = histcounts(lickTime(licksInbout), sort([boutIntvStart; boutIntvEnd]));
%     incintv = N(1:2:end) > minBoutLicks;

    % use 1 ms resolution
    % i don't need the full epoch timeseries because this eventually
    % becomes intervals again after being evaluated for a condition... 
    time = [boutIntvStart(1)-.001:.001:boutIntvEnd(end)+.001]';
    % this is faster but less flexible than list2vec
    lickBout = isIncluded(time, [boutIntvStart boutIntvEnd]);
    
    % create column timeFromLick that specifies time from closest lick
    [~, timeFromLick] = knnsearch(lickTime, time);
    
    lickBoutXPvec = isIncluded(time, [lickBoutXP lickBoutXP]);
    if output_intervals
        out{day}{epoch}.boutTimes = [boutIntvStart boutIntvEnd];
        out{day}{epoch}.lickBoutXP = lickBoutXP;
    else
        
        out{day}{epoch}.time =  time; % time index
        out{day}{epoch}.timeFromLick = timeFromLick;
        out{day}{epoch}.lickBout = lickBout; % binary, active lick bout or not
        out{day}{epoch}.lickBoutXP = lickBoutXPvec; % binary, XP within lick bout or not
    end
%     out{day}{epoch}.velocity = pos{day}{epoch}.data(:,9); % speed of animal
end











