

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
%subjects = {'S48','S55','S58','S59','S60','S61'};
subjects = {'S48'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for s = 1:numel(subjects)
  clear('spike_config');
  clear('spike');
  this_subject = subjects{s};

  try
    load(sprintf('%s/%s/%s_session.mat',path_prefix,this_subject,this_subject));
    assert(exist('session') == 1);
  catch
    error('could not load session struct array for %s',this_subject);
  end

  spike_config_filename = sprintf('%s/%s/spike/%s_spike_config.mat', ...
      path_prefix,this_subject,this_subject);
  try
    load(spike_config_filename);
    assert(exist('spike_config') == 1);
  catch
    error('could not load spike_config from %s',spike_config_filename);
  end

  for t = 1:length(spike_config)
    if ~strcmp(spike_config(t).subject,this_subject)
      % Skip if the subject doesn't match
      warning('spike_config(%d) does not match subject %s',t,this_subject);
      continue;
    end
    clear('clust');
    this_day = spike_config(t).day;

    % ugly hack
    if (this_day ~= 1)
      continue;
    end

    this_tetrode = spike_config(t).tetrode;
    this_depth = spike_config(t).depth;
    % Check that spike waveforms file exists, and skip if it doesn't
    spike_filename = sprintf( ...
        '%s/%s/spike/tetrode%02d/%s_day%d_tetrode%02d(%03d)_spike.mat', ...
        path_prefix,this_subject,this_tetrode,this_subject,this_day, ...
        this_tetrode,this_depth);
    if (exist(spike_filename) ~= 2)
      warning('expected spike waveforms file %s is not present; skipping', ...
          spike_filename);
      continue
    end
    % compute cluster parameters (this takes time!)
    try
      disp(['computing cluster parameters for ' spike_filename]);
      clust = make_matclust_struct(spike_filename, ...
          session(find([session(:).day] == this_day)));
    catch
      error('make_matclust_struct(%s,session) failed',spike_filename);
    end
    clust_filename = sprintf('%s_matclust.mat', ...
         regexp(spike_filename,'.+(?=_spike.mat$)','match','once'));
    try
      save(clust_filename,'-struct','clust');
    catch
      error('could not save matclust data as %s',clust_filename);
    end
  end
end
clear;

