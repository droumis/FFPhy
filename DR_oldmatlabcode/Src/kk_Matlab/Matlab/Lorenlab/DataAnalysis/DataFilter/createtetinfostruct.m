function createtetinfostruct(animdirect,fileprefix,append,varargin)
% createtetinfostruct(animdirect,fileprefix)
% createtetinfostruct(animdirect,fileprefix,append)
%
% This function creates a tetinfo file in the animal's data directory.
% For each tetrode, the depth and number of cells is saved.  If a tetinfo file
% exists and new data is being added, set APPEND to 1.

daydirect = '';
daynum = 0;
[otherArgs] = procOptions(varargin);

tetinfo = [];
if (nargin < 3)
   append = 0;
end
if (append)
   try
      load(fullfile(animdirect,[fileprefix,'tetinfo']));
   end
end


%create tetinfo for all tetrodes
if isempty(daydirect)
  eegfiles = dir(fullfile(animdirect, 'EEG', [fileprefix, 'eeg*']));
  if isempty(eegfiles)
     error('No EEG files found.');
  end

  for i = 1:length(eegfiles)
     load(fullfile(animdirect, 'EEG', eegfiles(i).name));
     depth =cellfetch(eeg, 'depth');
     try
        tetinfo{depth.index(1)}{depth.index(2)}{depth.index(3)}.depth = depth.values{1};
  %       tetinfo{depth.index(1)}{depth.index(2)}{depth.index(3)}.depth = depth.values;
        tetinfo{depth.index(1)}{depth.index(2)}{depth.index(3)}.numcells = 0;
        tetinfo{depth.index(1)}{depth.index(2)}{depth.index(3)}.tetrode = depth.index(3);
     catch
        warning([fileprefix, 'eeg0', num2str(depth.index(1)), '-', num2str(depth.index(2)), '-', num2str(depth.index(3)), '.mat',  ' may be empty'])
     end
  end
else
  eegfiles = dir(fullfile(daydirect,'*.eeg'));
  if isempty(eegfiles)
    error('No EEG files found.');
  end
  for i = 1:length(eegfiles)
    [tet,depth] = strread(eegfiles(i).name,'%d-%d.eeg');
    tetinfo{daynum}{tet}.depth = depth;
    tetinfo{daynum}{tet}.numcells = 0;
    tetinfo{daynum}{tet}.tetrode = tet;
  end
end

% add number of cells on each tetrode
spikefiles = dir(fullfile(animdirect, 'EEG', [fileprefix, 'spikes*']));

for i = 1:length(spikefiles)

   load(fullfile(animdirect, spikefiles(i).name));
   timerange = cellfetch(spikes, 'timerange');

   for j = 1:size(timerange.index,1)
      numcells = length(spikes{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)});
      tetinfo{timerange.index(j,1)}{timerange.index(j,2)}{timerange.index(j,3)}.numcells =numcells;
   end

end

save(fullfile(animdirect,[fileprefix,'tetinfo']),'tetinfo');
