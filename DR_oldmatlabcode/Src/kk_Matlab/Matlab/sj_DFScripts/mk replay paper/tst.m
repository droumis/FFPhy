cq = []
% Go through each set of cells from training filter and find each of the epochs
% where that cell has a set of cluster quality measures.
for a = 1:length(trainingfilter)
    for env = 1:length(trainingfilter(a).output)
	for epoch = 1:length(trainingfilter(a).output{env})
	    ind = trainingfilter(a).output{env}(epoch).index;
	    % find these cells in each epoch of qf
	    nqfe = length(qf(a).output)
	    q = ones(size(ind,1), nqfe) * -1;
	    for i = 1:length(qf(a).output)
		tmpi = rowfind(ind(:, [1 3 4]), qf(a).output{i}(:,[1 3 4]));
		found = find(tmpi);
		q(found,i) = qf(a).output{i}(tmpi(found),6);
	    end
	    cq = [cq; q];
	end
    end
end


