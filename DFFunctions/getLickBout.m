function out = getLickBout(animaldir, animal, epochs, varargin)

lickGap = .5; % upper inter bout-lick interval threshold in seconds 
boutNum = 10; %filter out bouts with less than boutNum licks
if ~isempty(varargin)
    assign(varargin{:})
end

% load data
days = unique(epochs(:,1));
lick = loaddatastruct(animaldir, animal, 'lick', days);
pos = loaddatastruct(animaldir, animal, 'pos', days);

for i = 1:size(epochs,1)
    day = epochs(i,1);
    epoch = epochs(i,2);
    
    % get lick bouts
    lickTimes = lick{day}{epoch}.starttime;
    boutIntvStart = lickTimes(find(diff([diff(lickTimes) < lickGap]) == 1)+1);
    boutIntvEnd = lickTimes(find(diff([diff(lickTimes) > lickGap]) == 1)+1);
    while boutIntvStart(1) > boutIntvEnd(1)
        boutIntvEnd(1) = [];
        if boutIntvEnd(end)<boutIntvStart(end)
            boutIntvStart(end) = [];
        end
    end
    % filter out bouts with less than boutNum licks
    licksInbout = logical(isExcluded(lickTimes, [boutIntvStart boutIntvEnd]));
    N = histcounts(lickTimes(licksInbout), sort([boutIntvStart; boutIntvEnd]));
    incintv = N(1:2:end) > boutNum;
    boutIntvStart = boutIntvStart(incintv);
    boutIntvEnd = boutIntvEnd(incintv);
    % downsample to the position times
    times = pos{day}{epoch}.data(:,1);
    lickBout = logical(isExcluded(times, [boutIntvStart boutIntvEnd]));
    
    % create column timeFromLick that specifies time from closest lick
    [~, timeFromLick] = knnsearch(lickTimes, times);
    
    out{day}{epoch}.time =  times; % time index
    out{day}{epoch}.timeFromLick = timeFromLick; 
    out{day}{epoch}.lickBout = lickBout; % binary, active lick bout or not
    out{day}{epoch}.velocity = pos{day}{epoch}.data(:,9); % speed of animal
end











