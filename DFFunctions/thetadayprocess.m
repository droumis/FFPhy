function thetadayprocess(animal,days, varargin)
%function thetadayprocess(rawdirectory,daydirectory,fileprefix,days, varargin)
%
%Applies a theta filter to all epochs for each day and saves the data in
%in the EEG subdirectory of the directoryname folder.  
%
%daydirectory - example '/data99/user/animaldatafolder/', a folder 
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
%        specifies, for each day, the tetrodes for which theta
%        extraction should be done
%
%     'downsample', factor
%        saves only every factor points (eg. 1:factor:end).
%                 Default 10
%
%     'f', matfilename
%        specifies the name of the mat file containing the
%        theta filter to use 
%        (default /usr/local/filtering/thetafilter.mat).  
%        Note that the filter must be called 'thetafilter'.

daytetlist = [];
f = '';
defaultfilter = 'thetadayprocess_filter.mat';
downsample = 1; % FILTER AUTOMATICALLY DOWNSAMPLES BY 10!!!
TSRATE = 10000;

[otherArgs] = procOptions(varargin);

% if the filter was not specified, load the default
if isempty(f) 
   eval(['load ', defaultfilter]);
else
   eval(['load ', f]);
end

animalinfo = animaldef(animal);
daydirectory = animalinfo{2};
fileprefix = animalinfo{3};

days = days(:)';

for day = days
  fprintf('[%d]',day);
  rawdirectory = animaldirdata(animal,day);
  % create the list of files for this day that we should filter
  tmpflist = dir(fullfile(rawdirectory,'*.eeg'));
  flist = cell(size(tmpflist));
  for i = 1:length(tmpflist)
    tetnum(i) = str2num(tmpflist(i).name(1:2));
    flist{i} = fullfile(rawdirectory,tmpflist(i).name);
  end
  if ~isempty(daytetlist)
    flist{~ismember(tetnum,daytetlist(find(daytetlist(:,1) == day),2))} = [];
    tetnum(~ismember(tetnum,daytetlist(find(daytetlist(:,1) == day),2))) = [];
    keyboard
  end

  if length(flist) == 0
    warning(sprintf('[thetadayprocess] No files found for day %d!',day));
    keyboard
  end

  timesdata = load(fullfile(rawdirectory,'times.mat'));
  NEpochs = size(timesdata.ranges,1)-1;

  % go through each file in flist and filter it
  for fnum = 1:length(flist) % (by tetrode)
    %load the eeg file
    for e = 1:NEpochs
      clear eeg_d;
      eeg_d.descript=flist{fnum};
      eeg_d.starttime=timesdata.ranges(e+1,1)/TSRATE;
      [eeg_d.data,eeg_t,eeg_d.samprate] = ...
        eeg_window(flist{fnum},timesdata.ranges(e+1,1)/TSRATE,diff(timesdata.ranges(e+1,:))/TSRATE);
      eeg_d.data = transpose(eeg_d.data);

      % filter it and save the result as single
      theta{day}{e}{tetnum(fnum)} = filtereeg2(eeg_d, ...
            thetafilter, 'single', 1); 
      thetafile = sprintf('%s/EEG/%stheta%02d-%d-%02d.mat', ...
               daydirectory, fileprefix, day, e, tetnum(fnum));
      % save the resulting file
      save(thetafile, 'theta');
      clear theta
    end % epochs
    fprintf('.');
  end % tetrode
  fprintf('\n');
end
