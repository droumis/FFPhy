function toppeak = dayprocess(daydirect, animdirect, fileprefix, daynum, varargin)

%Run this program from the folder containing all of the animal's day folders.  It will process the matclust files for the desired
%day, and put the files in the ANIMDIRECT directory.
%
%dayprocess(daydirect, animdirect, fileprefix, daynum, varargin)
%
%DAYDIRECT -- folder name where the day's data is stored
%ANIMDIRECT -- the path to where the animal's processed data will be stored -- example '/data99/student/Mil'
%FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
%DAYNUM -- the day number for the experiment (starting with 1)
%OPTIONS -- optional input in the form of: 'option', value, 'option', value
%           optional arguments are:
%           DIODEPOS -- 0 uses back diode for pos, 1 uses front diode for
%               pos, and values in between use the proportionate distance
%               between the two diodes. (default 0.5, ie, average of diodes)
%           CMPERPIX -- size of each pixel in centimeters (default 1)
%           POSFILT -- filter to use for pos smoothing when computing
%               velocity (default gaussian(30*0.5, 60))
%           PROCESSEEG -- 1 to process EEG, 0 to skip EEG processing (default 1)
%           SYSTEM -- 1 for old rig, 2 for nspike (default 2)
%           THETATETRODE -- -1 to use whichever tetrode had the largest averege theta, otherwise put the tetrode
%                           number to use as the theta reference for phase calculation (default -1)
%           VARPREFIX -- the first three letters of the animals name, which is attached to all variable names.
%              I recommend putting '' for this, which will not include any name in the variables. (Default '')
%           REVERSEX -- mirror image pos info along x axis (default 0)
%           REVERSEY -- mirror image pos info along y axis (default 0)

% processdio = 0;
% processstimdio = 1;
% processspikes = 1;
% processpos = 0;
% processrawpos = 1;

% set default values for the optional arguments
diodepos = 0.5;
cmperpix = 1;
processeeg = 1;
system = 2;
posfilt = gaussian(30*0.5, 60);
% thetatetrode = -1;
varprefix = '';
reversex = 0;
reversey = 0;
processdio = 0;
processstimdio = 0;
processspikes = 0;
processpos = 0;
processrawpos = 0;
processtask = 0;
processtetinfo = 0;

[otherOptions] = procOptions(varargin);

if (exist(fullfile(pwd,daydirect),'dir') == 7) % is day direct here now?
  daydirect = fullfile(pwd,daydirect);
end

[daypath,dayname] = fileparts(daydirect);

currdir = pwd;

[success,message] = mkdir(animdirect,'EEG');
if ~success
  error('Unable to make EEG directory. %s',message);
end

NUMTETRODES = 30;
load thetafilter;

currdir = pwd;
times = gettimes(fullfile(daydirect,'times.mat'));
runepochs = [];
eeg = [];

if (processpos)
  for epoch = 1:length(times)
    if ~isempty(times(epoch).starttime)
      starttime = times(epoch).starttime;
      endtime = times(epoch).endtime;
      epochname = times(epoch).name;
      % read in the position file
      cd(daydirect)
      posfiles = dir('*.p');


      % read in the raw diode information (front and back positions)
      % this looks through all the .p files and picks out the times that are
      % within the start and end times
      rawpos{daynum}{epoch} = readrawpos(posfiles,starttime,endtime);
      cd(currdir);

if 0
      % interpolate over the raw positions to get location and direction
      % at each time
      if ~isempty(rawpos{daynum}{epoch}.data)
        pos{daynum}{epoch} = posinterp(rawpos{daynum}{epoch}, 'diode', diodepos,...
          'maxv', 300, 'maxdevpoints', 30, 'reversex', reversex, 'reversey', reversey);

        % multiply pos data by cmperpix to convert from pixel units to cm
        pos{daynum}{epoch}.cmperpixel = cmperpix;
        pos{daynum}{epoch}.data(:,2:3) = pos{daynum}{epoch}.data(:,2:3)*cmperpix;

        % now run functions to add additional behavior parameters to pos struct
        pos{daynum}{epoch} = addvelocity(pos{daynum}{epoch}, posfilt);
      else
        pos{daynum}{epoch}.data = [];
        pos{daynum}{epoch}.cmperpixel = cmperpix;
      end

      if (strncmp(epochname, 'run', 3))
        runepochs = [runepochs epoch];
        % assign the environment and task variables
      end
end

    end
  end
  pos{daynum} = estimate_position(rawpos{daynum},'centimeters_per_pixel',cmperpix);

  eval([varprefix,'rawpos = rawpos;']);
  save(fullfile(animdirect,sprintf('%srawpos%02d.mat',fileprefix,daynum)),[varprefix,'rawpos']);
  eval([varprefix,'pos = pos;']);
  save(fullfile(animdirect,sprintf('%spos%02d.mat',fileprefix,daynum)),[varprefix,'pos']);
elseif (processrawpos)
  for epoch = 1:length(times)
    if ~isempty(times(epoch).starttime)
      epochname = times(epoch).name;
      % read in the position file
      posfiles = dir(fullfile(daydirect,'*rawpos*.mat'));

      % read in the raw diode information (front and back positions)
      % this looks through all the .p files and picks out the times that are
      % within the start and end times
      cd(daydirect)
      rawpos{daynum}{epoch} = readrawpos_steve(posfiles,times(epoch).range);
      cd(currdir)

if 0
      % interpolate over the raw positions to get location and direction
      % at each time
      if ~isempty(rawpos{daynum}{epoch}.data)
        pos{daynum}{epoch} = posinterp(rawpos{daynum}{epoch}, 'diode', diodepos,...
          'maxv', 300, 'maxdevpoints', 30, 'reversex', reversex, 'reversey', reversey);

        % multiply pos data by cmperpix to convert from pixel units to cm
        pos{daynum}{epoch}.cmperpixel = cmperpix;
        pos{daynum}{epoch}.data(:,2:3) = pos{daynum}{epoch}.data(:,2:3)*cmperpix;

        % now run functions to add additional behavior parameters to pos struct
        pos{daynum}{epoch} = addvelocity(pos{daynum}{epoch}, posfilt);
      else
        pos{daynum}{epoch}.data = [];
        pos{daynum}{epoch}.cmperpixel = cmperpix;
      end

      if (strncmp(epochname, 'run', 3))
        runepochs = [runepochs epoch];
        % assign the environment and task variables
      end
end
    end
  end
  pos{daynum} = estimate_position(rawpos{daynum},'centimeters_per_pixel',cmperpix);


  eval([varprefix,'rawpos = rawpos;']);
  save(fullfile(animdirect,sprintf('%srawpos%02d.mat',fileprefix,daynum)),[varprefix,'rawpos']);
  eval([varprefix,'pos = pos;']);
  save(fullfile(animdirect,sprintf('%spos%02d.mat',fileprefix,daynum)),[varprefix,'pos']);
else
  pos = loaddatastruct(animdirect,fileprefix,'pos',daynum);
end


if (processeeg) %this part filters the eeg data and saves the theta files (which takes some time)
  % get a list of all the eeg files in the directory
  eegfiles = dir(fullfile(daydirect,'*.eeg'));
  if ~isempty(eegfiles)
    eegfileinfo = regexp({eegfiles.name},'(?<tetrode>\d\d)-(?<depth>\d*)\.eeg','names');
    eegfileinfo = cell2mat(eegfileinfo);
    depth = str2num(cat(1,eegfileinfo.depth));
    eegtet = str2num(cat(1,eegfileinfo.tetrode));
    unique_tetrodes = unique(eegtet);
    fprintf('Processing %d eeg files.',length(eegfiles));
    for i = 1:length(unique_tetrodes)
      fprintf('.');
      t = unique_tetrodes(i);
      inds = find(eegtet==t);
      eegfilename = eegfiles(i).name;

      % now get the eeg data for each epoch
      for epoch = 1:length(times)
        if ~isempty(times(epoch).starttime)
          starttime = times(epoch).starttime;
          endtime = times(epoch).endtime;
          epochname = times(epoch).name;
          if (system == 1)
            eegstruct = readeeg(fullfile(daydirect,eegfiles(inds(1)).name), starttime, endtime, eegtet/NUMTETRODES);
          elseif (system == 2)
            for j = 1:length(inds)
              eegstruct = nreadeeg(fullfile(daydirect,eegfiles(inds(j)).name), starttime, endtime);
              if ~isempty(eegstruct)
                break;
              end
            end
          end
          eegstruct.depth = depth(inds(1));
          eval(sprintf('%seeg{%d}{%d}{%d} = eegstruct;', varprefix,daynum,epoch,eegtet(inds(1))));
          tmpfilename = fullfile(animdirect,'EEG',sprintf('%seeg%02d-%d-%02d.mat',fileprefix,daynum,epoch,eegtet(inds(1))));
          save(tmpfilename,[varprefix,'eeg']);
          eval(['clear ',varprefix,'eeg']);
        end
      end
    end
    fprintf('\n');
  end
end



if (processspikes)

  spikes = [];
  multi = [];

  tetfolders = dir(fullfile(daydirect,'*-*'));
  tetfolders(~[tetfolders.isdir]) = [];

  tetfolders = {tetfolders.name};
  tetfolderinfo = regexp(tetfolders,'(?<tetrode>\d\d)-(?<depth>\d\d\d)','names');

  for tet = 1:length(tetfolders)
  if ~isempty(tetfolderinfo{tet})
    depth = str2num(tetfolderinfo{tet}.depth);
    tetnum = str2num(tetfolderinfo{tet}.tetrode);

    if exist(fullfile(daydirect,tetfolders{tet},sprintf('%s-%02d_params.mat',dayname,tetnum)))
      rawdata = load(fullfile(daydirect,tetfolders{tet},sprintf('%s-%02d_params.mat',dayname,tetnum)));
      if isempty(rawdata)
        error(sprintf('Raw param file not found for tetrode %d.',tetnum));
      end
      for epoch = 1:length(times)
        if ~isempty(times(epoch).starttime)
          starttime = timetrans(times(epoch).starttime,10000,2);
          endtime = timetrans(times(epoch).endtime,10000,2);
          epochname = times(epoch).name;
          multi{daynum}{epoch}{tetnum} = rawdata.filedata.params((find((rawdata.filedata.params(:,1) >= starttime) & (rawdata.filedata.params(:,1) <= endtime))),1);
        end
      end
    else
      multi{daynum}{epoch}{tetnum} = [];
    end
    

    matclustfile = dir(fullfile(daydirect,tetfolders{tet},'matclust_*.mat'));

    if ~isempty(matclustfile)
      for epoch = 1:length(times)
        if ~isempty(times(epoch).starttime)
          starttime = times(epoch).starttime;
          endtime = times(epoch).endtime;
          epochname = times(epoch).name;
          clear clustattrib clustdata
          load(fullfile(daydirect,tetfolders{tet},matclustfile(1).name));
          % generate the list of times corresponding to the peaks of the theta rhythm
          starttime2 = timetrans({starttime},10000,2);
          endtime2 = timetrans({endtime},10000,2);
          if ~isempty(clustattrib.clusters)
            for clustnum = 1:length(clustattrib.clusters)
              if (~isempty(clustattrib.clusters{clustnum}))
                %make sure that the cluster was defined for the current time epoch.  If not, make the cell empty.
                if (is_cluster_defined_in_epoch(clustnum,epoch+1))
                  timepoints = clustdata.params(clustattrib.clusters{clustnum}.index,1);
                  if (size(clustdata.params,2) > 4)
                    amps = clustdata.params(clustattrib.clusters{clustnum}.index,2:5);
                  else
                    amps = nan(size(clustdata.params,1),1);
                  end
                  if (size(clustdata.params,2) > 5)
                    spikewidth = mean(clustdata.params(clustattrib.clusters{clustnum}.index,6));
                  else
                    spikewidth = nan(size(clustdata.params,1),1);
                  end

                  timepoints = timepoints(find((timepoints >= starttime2) & (timepoints <= endtime2)));
                  amps = amps(find((timepoints >= starttime2) & (timepoints <= endtime2)),:);
                  ampvar = var(amps);
                  [trash, maxvar] = max(ampvar);
                  amps = amps(:,maxvar);
                  timepoints = timepoints(:);
                  if isfield(pos{daynum}{epoch},'data') && ~isempty(pos{daynum}{epoch}.data)
                    %if the epoch was defined for the cluster, but no valid points were inside the boxes, make the data field empty
                    spikes{daynum}{epoch}{tetnum}{clustnum}.data = [];
                    spikes{daynum}{epoch}{tetnum}{clustnum}.descript = 'spike data';
                    spikes{daynum}{epoch}{tetnum}{clustnum}.fields = 'time x y dir thetaphase amplitude(highest variance channel) posindex';
                    spikes{daynum}{epoch}{tetnum}{clustnum}.depth = depth;
                    spikes{daynum}{epoch}{tetnum}{clustnum}.spikewidth = spikewidth;
                    spikes{daynum}{epoch}{tetnum}{clustnum}.timerange = [starttime2 endtime2];

                    [spikepos, posindex] = lookuptime3(timepoints/10000, pos{daynum}{epoch}.data(:,1),pos{daynum}{epoch}.data(:,2:4)');
                    findgoodpoints = find((spikepos(1,:)~=0)&(spikepos(2,:)~=0));
                    spikepos = spikepos(:,findgoodpoints);
                    posindex = posindex(findgoodpoints);
                    timepoints = timepoints(findgoodpoints);
                    amps = amps(findgoodpoints);

                    if ~isempty(timepoints)
                      spikes{daynum}{epoch}{tetnum}{clustnumdayprocess('/home/kkay/Documents/qayyim07282010','/home/kkay/Documents/qay','qay',1,'proc',0,'processeeg',1,'cmperpix',1)}.data = timepoints/10000;
                      spikes{daynum}{epoch}{tetnum}{clustnum}.data(:,2:4) = spikepos';
                      spikes{daynum}{epoch}{tetnum}{clustnum}.data(:,6) = amps;
                      spikes{daynum}{epoch}{tetnum}{clustnum}.data(:,7) = posindex;
                    end
                  else %no valid pos info, so just save the spikes
                    spikes{daynum}{epoch}{tetnum}{clustnum}.data = [timepoints/10000 amps];
                    spikes{daynum}{epoch}{tetnum}{clustnum}.descript = 'spike data';
                    spikes{daynum}{epoch}{tetnum}{clustnum}.fields = 'time amplitude(highest variance channel)';
                    spikes{daynum}{epoch}{tetnum}{clustnum}.depth = depth;
                    spikes{daynum}{epoch}{tetnum}{clustnum}.spikewidth = spikewidth;
                    spikes{daynum}{epoch}{tetnum}{clustnum}.timerange = [starttime2 endtime2];
                  end
                else
                  spikes{daynum}{epoch}{tetnum}{clustnum} = [];
                end
              end
            end
          end
        end
      end
    end
  end
  end
dayprocess('/home/kkay/Documents/qayyim07282010','/home/kkay/Documents/qay','qay',1,'proc',0,'processeeg',1,'cmperpix',1)
  eval([varprefix,'spikes = spikes;']);
  save(fullfile(animdirect,sprintf('%sspikes%02d.mat',fileprefix,daynum)),[varprefix,'spikes']);
  eval([varprefix,'multi = multi;']);
  save(fullfile(animdirect,sprintf('%smulti%02d.mat',fileprefix,daynum)),[varprefix,'multi']);

end


if (processdio)
  creatediostruct(daydirect,animdirect,fileprefix,daynum);
elseif (processstimdio)
  stimdiofiles = dir(fullfile(daydirect,'*stimdio*.mat'));
  if length(stimdiofiles) > 1
    error('[dayprocess: processstimdio] Multiple stimdio mat files found.');
  end
  tmp_stimdio = load(fullfile(daydirect,stimdiofiles.name),'stimdio');
  stimdio{daynum} = tmp_stimdio.stimdio;

  eval([varprefix,'stimdio = stimdio;']);
  save(fullfile(animdirect,sprintf('%sstimdio%02d.mat',fileprefix,daynum)),[varprefix,'stimdio']);
end

if (processtask)
  taskfiles = dir(fullfile(daydirect,'*task*.mat'));
  if length(taskfiles) == 1
    tmp_task = load(fullfile(daydirect,taskfiles.name),'task');
    task{daynum} = tmp_task.task;
    eval([varprefix,'task = task;']);
    save(fullfile(animdirect,sprintf('%stask%02d.mat',fileprefix,daynum)),[varprefix,'task']);
  end
end

if (processtetinfo)
  tetfiles = dir(fullfile(daydirect,'*tetinfo*.mat'));
  if length(tetfiles) == 1
    tmp_tet = load(fullfile(daydirect,tetfiles.name),'tetinfo');
    % tetinfo is a single file for all days and epochs
    % check to see if it already exists
    fname = fullfile(animdirect,sprintf('%stetinfo.mat',fileprefix));
    if exist(fname,'file')
      load(fname);
    end
    tetinfo{daynum} = tmp_tet.tetinfo;
    eval([varprefix,'tetinfo = tetinfo;']);
    save(fname,[varprefix,'tetinfo']);
  end
end


cd(currdir);
