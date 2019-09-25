function out = calcmuafiringrate(index, excludetimes, spikes, pos, varargin)
% out = calcfiringrate(index, excludetimes, spikes, options)
% Calculates an estimate of the population firing rate at all valid times.
%
% Options:
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%

smooth = [];
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'smooth'
                smooth = varargin{option+1};
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
timestep = median(diff(p.data(:,1)));
mua = zeros(size(p.data(:,1)));
time = p.data(:,1);


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

if ~isempty(excludetimes)
    goodtimes = logical(~isExcluded(p.data(:,1), excludetimes));
    time = time(goodtimes);
    mua = mua(goodtimes);
end

if ~isempty(smooth)
    g = gaussian(1,10);
    tmp = mua./timestep;
    out.mua = smoothvect(tmp,g); 
else
    out.mua = mua./timestep;
end
out.time = time;
out.numcells = numcells;
end