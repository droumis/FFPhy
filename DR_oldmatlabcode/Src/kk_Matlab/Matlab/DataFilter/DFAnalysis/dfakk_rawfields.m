function out = dfakk_rawfields(index, excludeperiods, spikes, pos)

% outputs non-excluded spikes and position data for each cell


spikedata = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
    includedspikes = spikedata(~isExcluded(spikedata(:,1), excludeperiods),1:3);
posdata = pos{index(1)}{index(2)}.data;
included_posindices = ~isExcluded(posdata(:,1),excludeperiods);
    included_posdata = posdata(included_posindices,1:3);

out.index = index;   
out.spikes = includedspikes;      % [time x y]
out.posdata = included_posdata;   % [time x y]

