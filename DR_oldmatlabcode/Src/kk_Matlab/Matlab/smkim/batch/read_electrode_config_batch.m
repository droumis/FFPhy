
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
continuous_hostnames = {'drizzle'};
spike_hostnames = {'rain'};
subjects = {'S48','S58','S59','S60','S61' };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)
  clear('continuous_config','spike_config','session','recording_sites');
  this_subject = subjects{s};

  % Read sessions meta-data and identify the days of recording
  session_filename = sprintf('%s/%s/%s_session.mat',path_prefix, ...
      this_subject,this_subject);
  try
    load(session_filename);
    assert(exist('session') == 1);
  catch
    error('could not load session file %s',session_filename);
  end
  if ~all(strcmp(this_subject,{session(:).subject}))
    error('not all elements of session match subject %s',this_subject);
  end  
  unique_days = unique([session(:).day]);

  % Construct a RECORDING_SITES data struct for this subject
  try
    recording_sites = read_recording_sites(sprintf( ...
      '%s/%s/%s_recording_sites.txt',path_prefix, ...
      this_subject,this_subject));
    assert(exist('recording_sites') == 1);
  catch
    error('error reading electrode positions inferred from histology');
  end
  if ~all(strcmp(this_subject,{recording_sites(:).subject}))
    error('not all elements of recording_sites match subject %s',this_subject);
  end  
  recording_sites_filename = sprintf('%s/%s/%s_recording_sites.mat', ...
      path_prefix,this_subject,this_subject);
  try
    disp(subjects{s});
    disp(unique({recording_sites(:).region}));
    save(recording_sites_filename,'recording_sites');
  catch
    error('could not save recording_sites as %s',recording_sites_filename);
  end

  %
  % For each day, call READ_CONTINUOUS_CONFIG and READ_SPIKE_CONFIG
  % with appropriate config filenames and recording_sites inputs
  for d = 1:numel(unique_days)
    this_day = unique_days(d);

    try
      tmp_continuous_config = read_continuous_config( ...
          cellfun(@(hostname) sprintf('%s/%s/%s_day%d.%s.config', ...
          path_prefix,this_subject,this_subject,this_day,hostname), ...
          continuous_hostnames,'UniformOutput',false), ...
          recording_sites(find([recording_sites(:).day] == this_day)));
    catch
      error('error reading continuous electrode config');
    end
    if (exist('continuous_config') == 1)
      continuous_config = [continuous_config; tmp_continuous_config];
    else
      continuous_config = tmp_continuous_config;
    end

    try
      tmp_spike_config = read_spike_config( ...
          cellfun(@(hostname) sprintf('%s/%s/%s_day%d.%s.config', ...
          path_prefix,this_subject,this_subject,this_day,hostname), ...
          spike_hostnames,'UniformOutput',false), ...
          recording_sites(find([recording_sites(:).day] == this_day)));
    catch
      error('error reading spike electrode config');
    end
    if (exist('spike_config') == 1)
      spike_config = [spike_config; tmp_spike_config];
    else
      spike_config = tmp_spike_config;
    end
  end

  % Save
  continuous_config_filename = sprintf('%s/%s/continuous/%s_continuous_config.mat', ...
      path_prefix,this_subject,this_subject);
  try
    save(continuous_config_filename,'continuous_config');
  catch
    error('could not save continuous_config as %s',continuous_config_filename);
  end
  spike_config_filename = sprintf('%s/%s/spike/%s_spike_config.mat', ...
      path_prefix,this_subject,this_subject);
  try
    save(spike_config_filename,'spike_config');
  catch
    error('could not save spike_config as %s',spike_config_filename);
  end
  %
 

end
clear;




