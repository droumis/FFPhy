function createcellinfostruct(animdirect,fileprefix,append)
% createcellinfostruct(animdirect,fileprefix)
% createcellinfostruct(animdirect,fileprefix,append)
%
% This function creates a cellinfo file in the animal's data directory.
% For each cell, the spikewidth, mean rate, and number of spikes is saved.  If a cellinfo file
% exists and new data is being added, set APPEND to 1.

cellinfo = [];
if (nargin < 3)
    append = 0;
end
if (append)
    try
        load(fullfile(animdirect,[fileprefix,'cellinfo']));
    end
end
spikefiles = dir(fullfile(animdirect,[fileprefix,'spikes*']));
for i = 1:length(spikefiles)
    load(fullfile(animdirect, spikefiles(i).name));
    timerange = cellfetch(spikes, 'timerange');
    for j = 1:size(timerange.index,1)
        cellinfo{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)}{timerange.index(j,4)}.spikewidth = timerange.values{j};
        try
            spikewidth = spikes{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)}{timerange.index(j,4)}.spikewidth;
        catch
            spikewidth = [];
        end
        epochtime = diff(spikes{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)}{timerange.index(j,4)}.timerange)/10000;
        numspikes = size(spikes{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)}{timerange.index(j,4)}.data,1);
        cellinfo{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)}{timerange.index(j,4)}.meanrate = numspikes/epochtime;
        cellinfo{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)}{timerange.index(j,4)}.numspikes = numspikes;
        cellinfo{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)}{timerange.index(j,4)}.spikewidth = spikewidth;
    end        
end

save(fullfile(animdirect,[fileprefix,'cellinfo']),'cellinfo');
