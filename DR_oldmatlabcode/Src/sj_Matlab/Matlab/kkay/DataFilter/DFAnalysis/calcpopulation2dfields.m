function out = calcpopulation2dfields(indices, excludeperiods, spikes, pos, binsize, peakthresh)
%trajdata = FILTERCALC2DFIELDS(index, excludeperiods, spikes, pos)
%trajdata = FILTERCALC2DFIELDS(index, excludeperiods, spikes, pos, binsize)
%
%Calculates the linear occupancy normalized firing rate for all the cells.
%
%spikes - the 'spikes' cell array for the day you are analyzing
%pos - the 'pos' cell array for the day you are analyzing
%indices - [day epoch tetrode cell]
%binsize- the length of each spatial bin (default 2cm)
%excludeperiods - [start end] times for each exlcude period
%
%The output is a structure. 


warning('OFF','MATLAB:divideByZero');
if (nargin < 5)
    binsize = 2;
end

i = [indices(1,1) indices(1,2)];


p = pos{i(1)}{i(2)};
% exclude the excluded positions
goodpos = ~isExcluded(p.data(:,1), excludeperiods);
p.data = p.data(goodpos,:);


out.ratemap = [];
out.minpos = [];
out.index = [];
validcount = 1;

%go through each cell and calculate the linearized rates
for cellcount = 1:size(indices,1)

    in = indices(cellcount,:);
    s = spikes{in(1)}{in(2)}{in(3)}{in(4)};
    
    if (~isempty(s.data))
	goodspikes = ~isExcluded(s.data(:,1), excludeperiods);
	s.data = s.data(goodspikes,:);
    end

    tmp = openfieldrate(s.data, p.data, binsize, 2);

    if max(tmp.spikerate >= peakthresh)
	out.ratemap{validcount} = tmp.smoothedspikerate;
	out.minpos = min(p.data(:,2:3));
	out.binsize = binsize;
        out.index = [out.index; in];
	validcount = validcount + 1;
    end
end

warning('ON','MATLAB:divideByZero');
