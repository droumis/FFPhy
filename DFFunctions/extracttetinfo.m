function extracttetinfo(animaldir, aname)
% extracttetinfo(animaldir, aname)
% load the cellinfo structure and create a tetinfo structure for this animal
cellinfo = loaddatastruct(animaldir, aname, 'cellinfo');
for d = 1:length(cellinfo)
    for e = 1:length(cellinfo{d})
	for t = 1:length(cellinfo{d}{e})
	    for c = 1:length(cellinfo{d}{e}{t})
		if (~isempty(cellinfo{d}{e}{t}{c}))
		    if (isfield(cellinfo{d}{e}{t}{c}, 'area'))
			tetinfo{d}{t}.area = cellinfo{d}{e}{t}{c}.area;
		    end
		end
	    end
	end
    end
end

save([animaldir, aname, 'tetinfo'], 'tetinfo');

