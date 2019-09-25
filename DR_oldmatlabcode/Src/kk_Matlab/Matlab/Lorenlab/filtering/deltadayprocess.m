function deltadayprocess(directoryname,fileprefix,days, varargin)
%THETADAYPROCESS(directoryname,fileprefix,days, options)
%
%Applies a delta filter to all epochs for each day and saves the data in
%in the EEG subdirectory of the directoryname folder.  
%
%directoryname - example '/data99/user/animaldatafolder/', a folder 
%            containing processed matlab data for the animal
%
%fileprefix -   animal specific prefix for each datafile (e.g. 'fre')
%
%days -         a vector of experiment day numbers 
%
%options -
%     'system', 1 or 2
%        specifies old (1) or new/nspike (2) rigs. Default 2.
%
%     'daytetlist', [day tet ; day tet ...]
%        specifies, for each day, the tetrodes for which delta
%        extraction should be done
%
%     'downsample', factor
%        saves only every factor points (eg. 1:factor:end).
%                 Default 10
%
%     'f', matfilename
%        specifies the name of the mat file containing the
%        delta filter to use 
%        (default /usr/local/filtering/deltafilter.mat).  
%        Note that the filter must be called 'deltafilter'.

daytetlist = [];
f = '';
defaultfilter = '';
system = 2;
downsample = 10;

[otherArgs] = procOptions(varargin);


minint = -32768;

% if the filter was not specified, load the default
if isempty(f) 
   eval(['load ', defaultfilter]);
else
   eval(['load ', f]);
end

days = days(:)';

for day = days
   % create the list of files for this day that we should filter
   if (isempty(daytetlist)) 
      tmpflist = dir(fullfile(directoryname,'EEG',sprintf('*eeg%02d-*.mat', day)));
      flist = cell(size(tmpflist));
      for i = 1:length(tmpflist)
         flist{i} = fullfile(directoryname,'EEG',tmpflist(i).name);
      end
   else
      % find the rows associated with this day
      flist = {};
      ind = 1;
      tet = daytetlist(find(daytetlist(:,1) == day),2);
      for t = tet;
         tmpflist = dir(fullfile(directoryname,'EEG', ...
           sprintf('*eeg%02d-*-%02d.mat',day,t)));
         nfiles = length(tmpflist); 
         for i = 1:length(tmpflist)
            flist{iind} = fullfile(directoryname,'EEG',tmpflist(i).name);
         end
      end
   end
   if length(flist) == 0
      warning(sprintf('[deltadayprocess] No files found for day %d!',day));
      keyboard
   end
   
   % go through each file in flist and filter it
   for fnum = 1:length(flist)
   % get the tetrode number and epoch
   % this is ugly, but it works
   dash = find(flist{fnum} == '-');
   epoch = str2num(flist{fnum}((dash(1)+1):(dash(2)-1)));
   tet = str2num(flist{fnum}((dash(2)+1):(dash(2)+3)));
      
   %load the eeg file
   load(flist{fnum});
   a = find(eeg{day}{epoch}{tet}.data < -30000);
   [lo,hi]= findcontiguous(a);  %find contiguous NaNs
   for i = 1:length(lo)
      if lo(i) > 1 & hi(i) < length(eeg{day}{epoch}{tet}.data)
      fill = linspace(eeg{day}{epoch}{tet}.data(lo(i)-1), ...
         eeg{day}{epoch}{tet}.data(hi(i)+1), hi(i)-lo(i)+1);
      eeg{day}{epoch}{tet}.data(lo(i):hi(i)) = fill;
      end
   end
      % filter it and save the result as %% int16, altered by kk 5.7.13
   delta{day}{epoch}{tet} = filtereeg2(eeg{day}{epoch}{tet}, ...
         deltafilter, 'int16', 1); 
      % NOTE THAT FILTER AUTOMATICALLY DOWNSAMPLES TO 150 Hz
      % downsample if requested 
   if (downsample ~= 1)
      delta{day}{epoch}{tet}.data = ...
      delta{day}{epoch}{tet}.data(1:downsample:end, :);
      delta{day}{epoch}{tet}.samprate = ...
         delta{day}{epoch}{tet}.samprate / downsample;
   end
   clear eegrec
   % replace the filtered invalid entries with the minimum int16 value of
   % -32768
   for i = 1:length(lo)
      if lo(i) > 1 & hi(i) < length(delta{day}{epoch}{tet}.data)
      delta{day}{epoch}{tet}.data(lo(i):hi(i)) = minint;
      end
   end


   % save the resulting file
   deltafile = sprintf('%s/EEG/%sdelta%02d-%d-%02d.mat', ...
            directoryname, fileprefix, day, epoch, tet);
   save(deltafile, 'delta');

   if 0
   % look for theta file to see if we should calculate T/D ratio
   fname = fullfile(directoryname,'EEG',sprintf('%stheta%02d-%d-%02d.mat',fileprefix,day,epoch,tet));
   if (exist(fname,'file'))
      % calculate T/D R
      load(fname);
      tdr{day}{epoch}{tet}.data = ...
         theta{day}{epoch}{tet}.data(:,2) ./ ...
         delta{day}{epoch}{tet}.data(:,2);

      tdrfile = sprintf('%s/EEG/%stdr%02d-%d-%02d.mat', ...
               directoryname, fileprefix, day, epoch, tet);
      save(tdrfile, 'tdr');
   end
   end
   
   clear delta
   end
end
