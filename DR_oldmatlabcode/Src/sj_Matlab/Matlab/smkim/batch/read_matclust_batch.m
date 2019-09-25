
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48', 'S58', 'S59', 'S60', 'S61'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)
  clear('spike_config');
  clear('unit');
  this_subject = subjects{s};
  spike_config_filename = sprintf('%s/%s/spike/%s_spike_config.mat', ...
      path_prefix,this_subject,this_subject);
  try
    load(spike_config_filename);
    assert(exist('spike_config') == 1);
  catch
    error('could not load spike_config from %s',spike_config_filename);
  end
  if ~all(strcmp(this_subject,{spike_config(:).subject}))
    error('spike_config does not match subject %s',this_subject);
  end
  unique_days = unique([spike_config(:).day]);
  for d = 1:numel(unique_days)
    this_day = unique_days(d);

    for t = 1:length(spike_config)
      if (spike_config(t).day ~= this_day)
        % skip to next iteration if the day doesn't match
        continue;
      end
      this_tetrode = spike_config(t).tetrode;
      this_depth = spike_config(t).depth;
      % Check that matclust file exists, and skip if it doesn't
      matclust_filename = sprintf( ...
          '%s/%s/spike/tetrode%02d/%s_day%d_tetrode%02d(%03d)_matclust.mat', ...
          path_prefix,this_subject,this_tetrode,this_subject,this_day, ...
          this_tetrode,this_depth);
      if (exist(matclust_filename) ~= 2)
        warning('expected matclust file %s is not present; skipping', ...
            matclust_filename);
        continue
      end
      % Read matclust data into struct array
      try
        tmp_unit = read_matclust(matclust_filename);
      catch
        error('read_matclust(%s) failed',matclust_filename);
      end
      if (exist('unit') ~= 1)
        unit = tmp_unit;
      else
        unit = [unit; tmp_unit];
      end
    end
    assert(all([unit(:).day] == this_day));
    % save combined UNIT struct array for this subject/day
    try
      save(sprintf('%s/%s/spike/%s_day%d_unit.mat',path_prefix, ...
          this_subject,this_subject,this_day),'unit');
    catch
      error('could not save single-unit spike data for subject %s day %d', ...
        this_subject,this_day);
    end
    clear('unit');
  end
end
clear;

