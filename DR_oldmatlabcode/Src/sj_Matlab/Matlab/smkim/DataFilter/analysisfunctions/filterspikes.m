function out = filterspikes(ind, excludetimes, spikes)

epochind = ind(:,1:2);
epochind = unique(epochind,'rows');
if (size(epochind,1) ~= 1)
    error('Indices must have only one unique day/epoch');
end

for i = 1:size(ind,1)
    if ~isempty(spikes{ind(i,1)}{ind(i,2)}{ind(i,3)}{ind(i,4)})
        if ~isempty(spikes{ind(i,1)}{ind(i,2)}{ind(i,3)}{ind(i,4)}.data)
            out = spikes{ind(i,1)}{ind(i,2)}{ind(i,3)}{ind(i,4)};
	    out.data = out.data(find(~isExcluded(out.data(:,1), ...
		    excludetimes)),:);
	end
    end
end
