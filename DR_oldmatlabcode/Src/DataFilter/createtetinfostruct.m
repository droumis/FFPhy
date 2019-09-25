function createtetinfostruct(animdirect,fileprefix)
% createtetinfostruct(animdirect,fileprefix)
% createtetinfostruct(animdirect,fileprefix,append)
%
% This function creates a tetinfo file in the animal's data directory.
% For each tetrode, the depth and number of cells is saved.  If a tetinfo file
% exists and new data is being added, set APPEND to 1.

if (animdirect(end) ~= '/')
    animdirect = [animdirect '/'];
end

tetinfo = [];
if (nargin < 3)
    append = 0;
end
if (append)
    try
        load([animdirect, fileprefix,'tetinfo']);
    end
end


%create tetinfo for all tetrodes

eegfiles = dir([animdirect, 'EEG/', fileprefix, 'eeg*']);
spikefiles = dir([animdirect, fileprefix, 'spikes*']);

 for i = 1:length(eegfiles)
    
     load([animdirect, 'EEG/', eegfiles(i).name]);
     depth =cellfetch(eeg, 'depth');
   
     try
       tetinfo{depth.index(1)}{depth.index(2)}{depth.index(3)}.depth = depth.values;
       tetinfo{depth.index(1)}{depth.index(2)}{depth.index(3)}.numcells = 0;
   catch
       warning([fileprefix, 'eeg0', num2str(depth.index(1)), '-', num2str(depth.index(2)), '-', num2str(depth.index(3)), '.mat',  ' may be empty'])
   end
   
 end
   
% add number of cells on each tetrode

for i = 1:length(spikefiles)
    
    load([animdirect, spikefiles(i).name]);
    timerange = cellfetch(spikes, 'timerange');

    for j = 1:size(timerange.index,1)
        numcells = length(spikes{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)});
        tetinfo{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)}.numcells =numcells;
    end
    
end

save([animdirect,fileprefix,'tetinfo'],'tetinfo');