%{
get the lick bout intervals
inputs:
epochs = n x 2 where n = [day epoch]

% what scripts created the lick data structure? i need to save the timeseries 
% dio times along with..
/home/droumis/Src/Matlab/filterframework_dr/Functions/get_licks.m

%}

function out = getLickBout(animal, epochs, varargin)

maxIntraBurstILI = .25; % max burst ili threshold in seconds
minBoutLicks = 3; %filter out bouts with less than boutNum licks
lick = [];
if ~isempty(varargin)
    assign(varargin{:})
end

% load data
days = unique(epochs(:,1));
if isempty(lick)
    andef = animaldef(animal);
    lick = loaddatastruct(andef{2}, animal, 'lick', days);
end
% pos = loaddatastruct(animaldir, animal, 'pos', days);
out = {};
for i = 1:size(epochs,1)
    day = epochs(i,1);
    epoch = epochs(i,2);
    
    % get lick bouts
    try
        lickTime = lick{day}{epoch}.eventtime;
    catch
        lickTime = lick{day}{epoch}.starttime; % legact name
    end
    ili = diff(lickTime);
    boutIntvStart = lickTime(find(diff([ili < maxIntraBurstILI]) == 1)+1);
    boutIntvEnd = lickTime(find(diff([ili > maxIntraBurstILI]) == 1)+1);
    while boutIntvStart(1) > boutIntvEnd(1)
        boutIntvEnd(1) = [];
        if boutIntvEnd(end)<boutIntvStart(end)
            boutIntvStart(end) = [];
        end
    end
    % filter out bouts with less than minBoutLicks
    try
        licksInbout = logical(isExcluded(lickTime, [boutIntvStart boutIntvEnd]));
    catch
        fprintf('error defining lick bouts for %d %d\n', day, epoch)
        continue
    end
    N = histcounts(lickTime(licksInbout), sort([boutIntvStart; boutIntvEnd]));
    incintv = N(1:2:end) > minBoutLicks;
    boutIntvStart = boutIntvStart(incintv);
    boutIntvEnd = boutIntvEnd(incintv);
    % use 1 ms resolution
    % i don't need the full epoch timeseries because this eventually
    % becomes intervals again after being evaluated for a condition... 
    time = [boutIntvStart(1)-.001:.001:boutIntvEnd(end)+.001]';
    % this is faster but less flexible than list2vec
    lickBout = logical(isExcluded(time, [boutIntvStart boutIntvEnd]));
    
    % create column timeFromLick that specifies time from closest lick
    [~, timeFromLick] = knnsearch(lickTime, time);
    
    out{day}{epoch}.time =  time; % time index
    out{day}{epoch}.timeFromLick = timeFromLick; 
    out{day}{epoch}.lickBout = lickBout; % binary, active lick bout or not

%     out{day}{epoch}.velocity = pos{day}{epoch}.data(:,9); % speed of animal
end











