function newdayprocess(animdir, animprefix, daynum, varargin)

% NOTE that if you used Loren's position reconstruction (.p files), use
% process_p
% if you used Steve's, then use process_pos


daydirect = pwd;          %%%%%%%%%%% change manually         
process_pos = 0;
process_pos_p = 0;
process_stim = 0;
process_stimeeg = 0;
process_longstimeeg = 1;     %% changed on 7.18.11 to extract windows

overwrite = 0;

FS = 30000;


% options for stimdio processing
%% note that these units are in seconds
longstimeeg_windowLength = 5;
longstimeeg_winOffset = 0;
longstimeeg_downsample = [15];        % 30 kHz to 2 kHz is a factor of 15
stimeeg_windowLength = 0.3;
stimeeg_winOffset = 0.1;

refelect = {};


[otherOptions] = procOptions(varargin);

if isempty(dir(fullfile(daydirect,'times.mat')))
  error('times.mat file not found in daydirect %s', daydirect);
end
% Load times file
times = gettimes(fullfile(daydirect,'times.mat'));
n_epochs = numel(times);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process Stimstructs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (process_stim)
  stimdio_filename = fullfile(animdir, sprintf('%sstimdio%02d.mat',animprefix, daynum));
  if (overwrite == 0) & ~isempty(dir(stimdio_filename))
    fprintf('Stimdio file %s exists, overwrite flag not set.', stimdio_filename);
  else
    %% Look for stimstructs and yank eeg windows
    stim_file = dir(fullfile(daydirect,'stimdio.mat'));
    if isempty(stim_file)
      error(sprintf('stimdio.mat file not found in daydirect (%s)\n',daydirect));
    end
    stim_file_contents = load(fullfile(daydirect,stim_file.name));
    stimdio{daynum} = stim_file_contents.stimdio;                         %% stimdio
    if numel(stimdio{daynum}) < n_epochs
      stimdio{daynum}{n_epochs} = [];                                     %% stimdio size initialization
    end

    if ~isempty(dir(stimdio_filename))
      fprintf('Overwrite flag set: overwriting %s\n', stimdio_filename);
    end
    %% Save stimdio files
    save(stimdio_filename,'stimdio');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process Stimeegs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (process_stimeeg)
  stimdio_filename = fullfile(animdir, sprintf('%sstimdio%02d.mat',animprefix, daynum));

  load(stimdio_filename);
  if ~exist('stimdio')
    error(sprintf('stimdio variable not found in %s\n', stimdio_filename));
  end

  stimeeg_filename = fullfile(animdir, sprintf('%sstimeeg%02d.mat',animprefix, daynum));   %% sets up where to save stimeeg
  
  if (overwrite == 0) & ~isempty(dir(stimeeg_filename))
    fprintf('Stimeeg file %s exists, overwrite flag not set.', stimeeg_filename);
  else
    i = 1;
    while isempty(stimdio{daynum}{i})
      i = i + 1;
    end
    n_pins = numel(stimdio{daynum}{i});

    %% Extract eeg windows around stimulation times
    t = [0:FS*stimeeg_windowLength-1]/FS - stimeeg_winOffset;
    cd(daydirect)
    d = dir('*.cont');
    if isempty(d)
        error('No CONT files found.');
    end
    cd(currdir);

       %% extracts data from stimdio >> stimeeg  , essentially replicating
       %% the struct structure
    for p = 1:n_pins
      for e = 1:n_epochs
        filename = d(i).name;
        if ~isempty(stimdio{daynum}{e}) && ~isempty(stimdio{daynum}{e})
          stimeeg{daynum}{e}.pulsetime = stimdio{daynum}{e}.pulsetimes(:,1);
          stimeeg{daynum}{e}.pulselength = stimdio{daynum}{e}.pulselength;
          stimeeg{daynum}{e}.data = nan(length(t),size(stimdio{daynum}{e}.pulsetimes,1),length(d));
          stimeeg{daynum}{e}.pin = stimdio{daynum}{e}.pin;
          stimeeg{daynum}{e}.type = stimdio{daynum}{e}.type;
          stimeeg{daynum}{e}.t = t;
        end
      end
    end

    for i = 1:length(d)
      for p = 1:n_pins
        for e = 1:n_epochs
          filename = d(i).name;
          if ~isempty(stimdio{daynum}{e}) && ~isempty(stimdio{daynum}{e})
            ts = stimdio{daynum}{e}.pulsetimes(:,1);                         %% timestamps of all stimulation pulses
            [stimeeg{daynum}{e}.data(:,:,i)] = cont_window_c(filename, ...   %% transfers data in the eeg window
              ts/10000 - stimeeg_winOffset, stimeeg_windowLength);
          end
        end
      end
    end

    
    %% if you want to UN-reference an electrode, specify the refelect above
    %% and this block of code will add the reference data back in
    if ~isempty(refelect)
      for i = 1:length(d)
        for p = 1:n_pins
          for e = 1:n_epochs
            filename = d(i).name;
            if ~isempty(stimdio{daynum}{e}) && ~isempty(stimdio{daynum}{e})
              for r = 1:length(refelect{i})
                stimeeg{daynum}{e}.data(:,:,i) = ...
                  stimeeg{daynum}{e}.data(:,:,i) + stimeeg{daynum}{e}.data(:,:,refelect{i}(r));  %% de-referencing
              end
            end
          end
        end
      end
    end

    %% Save stimeeg files
    if ~isempty(dir(stimeeg_filename))
      fprintf('Overwrite flag set: overwriting %s\n', stimeeg_filename);
    end
    save(stimeeg_filename,'stimeeg');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process Long Stimeegs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (process_longstimeeg)
  stimdio_filename = fullfile(animdir, sprintf('%sstimdio%02d.mat',animprefix, daynum));

  load(stimdio_filename);
  if ~exist('stimdio')
    error(sprintf('stimdio variable not found in %s\n', stimdio_filename));
  end

  longstimeeg_filename = fullfile(animdir, sprintf('%slongstimeeg%02d.mat',animprefix, daynum));  %% sets up longstimeeg
  
  if (overwrite == 0) & ~isempty(dir(longstimeeg_filename))
    fprintf('Stimeeg file %s exists, overwrite flag not set.', stimeeg_filename);
  else
    i = 1;
    while isempty(stimdio{daynum}{i})
      i = i + 1;
    end
    n_pins = numel(stimdio{daynum}{i});

    %% Extract eeg windows around stimulation times
    t = [0:FS*longstimeeg_windowLength-1]/FS - longstimeeg_winOffset;   %% vector of times for a given window
    for s = 1:length(longstimeeg_downsample)
      t = decimate(t, longstimeeg_downsample(s));                       %% downsampled times
    end
    currdir = pwd;
    cd(animdir);                %% changed from daydirect to animdirect
    d = dir('*.cont');
    if isempty(d)                                       %% checks for .cont file
        error('No CONT files found.');
    end
    cd(currdir);

    
        %% as in stimeeg, this extracts the basic fields
   % for p = 1:n_pins          %% removed all code relating to pin struct
      for e = 1:n_epochs
        filename = d(i).name;
        if ~isempty(stimdio{daynum}{e})
          longstimeeg{daynum}{e}.pulsetime = stimdio{daynum}{e}.pulsetimes(:,1);
          longstimeeg{daynum}{e}.pulselength = stimdio{daynum}{e}.pulselength;
          longstimeeg{daynum}{e}.data = nan(length(t),size(stimdio{daynum}{e}.pulsetimes,1),length(d));
          longstimeeg{daynum}{e}.pin = stimdio{daynum}{e}.pin;
          longstimeeg{daynum}{e}.type = stimdio{daynum}{e}.type;
          longstimeeg{daynum}{e}.t = t;
        end
      end
 %   end

    % extract eeg windows
    for i = 1:length(d)
      for p = 1:n_pins
        for e = 1:n_epochs
          filename = d(i).name;
          if ~isempty(stimdio{daynum}{e}) && ~isempty(stimdio{daynum}{e})
            ts = stimdio{daynum}{e}.pulsetimes(:,1);                             %% timestamps of stimulation pulses
            ex = cont_window_c(filename, ts/10000 - longstimeeg_winOffset, longstimeeg_windowLength);
            for s = 1:length(longstimeeg_downsample)
              ex = matrix_decimate(ex, longstimeeg_downsample(s));                  %% downsampling of eeg data
            end
            longstimeeg{daynum}{e}.data(:,:,i) = ex;
          end
        end
      end
    end

    %% de-referencing, if specified
    if ~isempty(refelect)
      for i = 1:length(d)
        for p = 1:n_pins
          for e = 1:n_epochs
            filename = d(i).name;
            if ~isempty(stimdio{daynum}{e}) && ~isempty(stimdio{daynum}{e})
              for r = 1:length(refelect{i})
                longstimeeg{daynum}{e}.data(:,:,i) = ...
                  longstimeeg{daynum}{e}.data(:,:,i) + longstimeeg{daynum}{e}.data(:,:,refelect{i}(r));
              end
            end
          end
        end
      end
    end

    %% Save longstimeeg files
    if ~isempty(dir(longstimeeg_filename))
      fprintf('Overwrite flag set: overwriting %s\n', longstimeeg_filename);
    end
    save(longstimeeg_filename,'longstimeeg');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process EEGs - convert to 1500 Hz sampling rate, save into EEG .mat files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process Spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process position data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (process_pos)
  fprintf('Processing rawpos position data...\n');
  %%% Look for rawpos files generated by Steve's code
  posfiles = dir(fullfile(daydirect,'rawpos*.mat')); % find rawpos files
  currdir = pwd;
  cd(daydirect)
  for epoch = 1:length(times)
    if ~isempty(times(epoch).starttime)
      rawpos{daynum}{epoch} = readrawpos_steve(posfiles,times(epoch).range);
    end
  end
  cd(currdir);

  %%% Look for .p files

  %%%%%%%%

  %% Save rawpos files
  save(fullfile(animdir, sprintf('%srawpos%02d.mat',animprefix, daynum)),'rawpos');

  %% Convert rawpos to pos using the Steve smoothing code
  pos{daynum} = estimate_position(rawpos{daynum},otherOptions{:});
  save(fullfile(animdir, sprintf('%spos%02d.mat',animprefix, daynum)),'pos');

  %% Smooth velocity with a gaussian kernel
  velocitydayprocess(animdir,animprefix,daynum);
end

if (process_pos_p)
  fprintf('Processing p-file position data...\n');
  %%% Look for .p files
  posfiles = dir(fullfile(daydirect,'*.p')); % find rawpos files
  currdir = pwd;
  cd(daydirect)
  for epoch = 1:length(times)
    if ~isempty(times(epoch).starttime)
      rawpos{daynum}{epoch} = readrawpos(posfiles,times(epoch).starttime,times(epoch).endtime);
    end
  end

  %%%%%%%%

  %% Save rawpos files
  save(fullfile(animdir, sprintf('%srawpos%02d.mat',animprefix, daynum)),'rawpos');

  %% Convert rawpos to pos using the Steve smoothing code
  pos{daynum} = estimate_position(rawpos{daynum},otherOptions{:});
  save(fullfile(animdir, sprintf('%spos%02d.mat',animprefix, daynum)),'pos');

  %% Smooth velocity with a gaussian kernel
  velocitydayprocess(animdir,animprefix,daynum);
end



function out = matrix_decimate(in, rate)
for i = 1:size(in,2)
  out(:,i) = decimate(in(:,i),rate);
end
