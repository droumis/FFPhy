function out = calcspeedfiringrate(index, excludetimes, spikes, pos, varargin)
% out = calcfiringrate(index, excludetimes, spikes, options)
% Calculates an estimate of the population firing rate at all valid times.
%
% Options:
%   'timebin', defines the windows to compute speed and firing rate.
%           Default = 0.5 seconds.
%

timebin = 0.5;

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'timebin'
                timebin = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

epochind = index(:,1:2);
epochind = unique(epochind,'rows');
if (size(epochind,1) ~= 1)
    error('Indices must have only one unique day/epoch');
end

p = pos{epochind(1)}{epochind(2)};
time = p.data(:,1);
mua = zeros(size(time));

% Determine how many spikes occur at all time for all cells
s = spikes{epochind(1)}{epochind(2)};
numcells = 0;
for i = 1:size(index,1)
    if ~isempty(s{index(i,3)}{index(i,4)})
        numcells = numcells+1;
        if ~isempty(s{index(i,3)}{index(i,4)}.data)
            mua(s{index(i,3)}{index(i,4)}.data(:,7)) = mua(s{index(i,3)}{index(i,4)}.data(:,7)) + 1;
        end
    end
end

% Apply excludetimes
if ~isempty(excludetimes)
    goodtimes = logical(~isExcluded(p.data(:,1), excludetimes));
    time = time(goodtimes);
    mua = mua(goodtimes);
end

% Determine speed and firing rate in each timebin
bin = lookup(time(1):timebin:time(end),time);
speed = nan(size(bin));
firing = nan(size(bin));

for s = 1:length(bin)
    if s < length(bin)
        speed(s) = mean(p.data(bin(s):bin(s+1),8));
        firing(s) = sum(mua(bin(s):bin(s+1)))./timebin;
    elseif s==length(bin)
        speed(s) = mean(p.data(bin(s):end,8));
        firing(s) = sum(mua(bin(s):end))./timebin;
    end    
end

out.firing = firing;
out.speed = speed;
out.numcells = numcells;

end