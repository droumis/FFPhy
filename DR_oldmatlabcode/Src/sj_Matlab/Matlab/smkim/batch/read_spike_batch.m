
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
%subjects = {'S48','S55','S58','S59','S60'};
subjects = {'S61'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first sweep through all folders to establish amplitude_cutoff for each .tt
% file, which we store in an args struct array

% args is a global struct array for storing arguments with which to call
% READ_SPIKE
args = struct('tt_filename',{},'spike_config',{},'amplitude_cutoff',{});

for s = 1:numel(subjects)
  clear('session');
  clear('spike_config');
  clear('spike');
  this_subject = subjects{s};
  try
    load(sprintf('%s/%s/%s_session.mat',path_prefix,this_subject,this_subject));
    assert(exist('session') == 1);
  catch
    error('could not load session info for subject %s',this_subject);
  end
  if ~all(strcmp(this_subject,{session(:).subject}))
    error('not all session info correspond to the same subject %s', ...
        this_subject);
  end
  try
    load(sprintf('%s/%s/spike/%s_spike_config.mat',path_prefix, ...
      this_subject,this_subject));
    assert(exist('spike_config') == 1);
  catch
    error('could not load spike_config for subject %s',this_subject);
  end
  % Pre-allocate args for this subject
  this_args = repmat( ...
      struct('tt_filename',{},'spike_config',{},'amplitude_cutoff',{}), ...
      size(spike_config));
  for t = 1:length(spike_config)
    % Check that the subject field of spike_config matches this_subject
    if ~strcmp(spike_config(t).subject,this_subject)
      error('subject field of spike_config does not match expected subject');
      continue;
    end
    this_day = spike_config(t).day;
    this_session = session(find([session(:).day] == this_day));
    this_tetrode = spike_config(t).tetrode;
    this_depth = spike_config(t).depth;
    tt_filename = sprintf('%s/%s/spike/tetrode%02d/%s_day%d_tetrode%02d(%03d).tt', ...
        path_prefix,this_subject,this_tetrode,this_subject,this_day,this_tetrode,this_depth);
    % Allow user to interactively set amplitude_cutoff
    try
      amplitude_cutoff = check_tt(tt_filename,this_session);
    catch
      error('check_tt(%s) failed',tt_filename);
    end
    this_args(t,1).tt_filename = tt_filename;
    this_args(t,1).spike_config = spike_config(t);
    this_args(t,1).amplitude_cutoff = amplitude_cutoff;
  end

  % Append this subject's args to the global args struct
  if isempty(args)
    args = this_args;
  else
    args = [args; this_args];
  end
end

% now convert the selected tt data to matlab structs
disp('now reading waveforms (this will take some time');
% then, sweep through again to run READ_SPIKE
for i = 1:numel(args)
  if isempty(args(i).tt_filename) || any(args(i).amplitude_cutoff == 0)
    continue; % skip if there are no waveforms to import
  end
  % import spike data from NSpike *.tt file
  try
    spike = read_spike(args(i).tt_filename,args(i).spike_config, ...
        args(i).amplitude_cutoff);
  catch
    error('failed to read spike waveforms from %s',args(i).tt_filename);
  end
  % skip if the spike data struct contains no events        
  if isempty(spike.timestamp)
    warning(['spike data struct will not be saved because it contains ' ...
        'no trigger events']);
    continue
  end
  [path, basename] = fileparts(args(i).tt_filename);
  spike_filename = sprintf('%s/%s_spike.mat',path,basename); ...
  % spike data struct gets saved in the individual tetrode folders
  try
    save(spike_filename,'spike');
  catch
    error('could not save spike waveforms as %s',spike_filename);
  end
  clear('spike');
end

clear;

