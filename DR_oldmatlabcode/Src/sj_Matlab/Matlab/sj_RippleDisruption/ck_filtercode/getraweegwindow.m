function [out] = getraweegwindow(index, times, excludeperiods, eegfiledir, varargin)
% function [out] = getraweegwindow(index, times, excludeperiods, eegfiledir, varargin)
%
%   index [day epoc tetrode]
%
%   out is windowed eeg organized as trial x time x tetrode

% assign the options

eegwindow = [-0.25 0.5];
[otherArgs] = procOptions(varargin);

out = [];
if isempty(times) | isempty(index)
  return;
end

eegfiles = dir(fullfile(eegfiledir,'*.eeg'));
for t = 1:length(eegfiles)
  [tet,depth] = strread(eegfiles(t).name,'%d-%d.eeg');
  tetrodedata(t).filename = eegfiles(t).name;
  tetrodedata(t).tetrode = tet;
end

for t = 1:size(index,1) % number of tetrodes
  goodPulses = ~isExcluded(times(:,1), excludeperiods);
  if size(times,2) > 1
    for tt = 1:size(times,2)
      goodPulses = goodPulses & ~isExcluded(times(:,1), excludeperiods);
    end
  end
  goodPulses = find(goodPulses);

  starttimes = times(goodPulses,1);

  starttimes = starttimes + eegwindow(1);
  
  eegfile = fullfile(eegfiledir,tetrodedata([tetrodedata(:).tetrode]==index(t,3)).filename);

  [out(:,:,t), actualtimes] = eeg_window(eegfile, starttimes, eegwindow(2)-eegwindow(1));
  % keyboard
end

