
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
%subjects = {'S48', 'S58', 'S59', 'S60', 'S61'};
subjects = {'S61'}

% load filters
load('theta_filter');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)
  clear('continuous_config');
  this_subject = subjects{s};
  try
    load(sprintf('%s/%s/%s_session.mat',path_prefix, ...
        this_subject,this_subject));
    assert(exist('session') == 1);
  catch
    error('could not load session info for %s',this_subject);
  end
  continuous_config_filename = sprintf( ...
      '%s/%s/continuous/%s_continuous_config.mat', ...
      path_prefix,this_subject,this_subject);
  try
    load(continuous_config_filename);
    assert(exist('continuous_config') == 1);
  catch
    error('could not load continuous_config from %s', ...
        continuous_config_filename);
  end
  if ~all(strcmp(this_subject,{continuous_config(:).subject}))
    error('continuous_config does not match subject %s',this_subject);
  end

  % iterate over days for this subject
  %unique_days = unique([continuous_config(:).day]);
  unique_days = 2;
  for d = 1:numel(unique_days)
    this_day = unique_days(d);
    clear('lfp');
    clear('emg');

    for t = 1:length(continuous_config)
      clear('continuous');
      % check all entries of continuous_config that match this_day
      if (continuous_config(t).day ~= this_day)
        continue;
      end
      this_electrode = continuous_config(t).electrode;
      this_channel = continuous_config(t).channel;
      this_depth = continuous_config(t).depth;
      % Check that file exists, and raise warning if it doesn't
      filename = sprintf( ...
          '%s/%s/continuous/%s_day%d_electrode%02d(%03d)_channel%d.eeg', ...
          path_prefix,this_subject,this_subject,this_day,this_electrode, ...
          this_depth,this_channel);
      if (exist(filename) ~= 2)
        warning('expected file %s is not present',filename);
        continue
      end
      % Read data into continuous data struct
      try
        continuous = read_continuous(filename,continuous_config(t), ...
            session([session(:).day] == this_day));
      catch
        error('read_continuous failed for %s',filename);
      end

      if strcmp('nuchal EMG',continuous_config(t).region)
        % append to end of emg
        if (exist('emg') == 1)
          emg = [emg; continuous];
        else
          emg = continuous;
        end
        continue;
      end
  
      % Determine the reference
      refelect = unique(vertcat(continuous(:).reference),'rows');
      if (size(refelect,1) ~= 1)
        error('references for %s are inconsistent',filename);
      end
      % skip if this continuous channel is already referenced to ground
      if isequal(refelect,[0 0]);
        % pass
      else
        % find the element of continuous_config that corresponds to the reference
        % channel
        u = find( ([continuous_config(:).day] == this_day) & ...
            ([continuous_config(:).electrode] == refelect(1)) & ...
            ([continuous_config(:).channel] == refelect(2)) );
        if (numel(u) > 1)
          error(['more than one element of continuous_config matches ' ...
              'reference channel (%s) for %s'],mat2str(refelect),filename);
        end
        if (numel(u) < 1)
          warning(['continuous_config specifies reference channel %s for %s ' ...
              'but this channel was not recorded!'], ...
              mat2str(refelect),filename);
          u = find( ([continuous_config(:).day] == this_day) & ...
            ([continuous_config(:).electrode] == refelect(1)) );
          if numel(u) == 1
            warning('instead will use reference electrode %d channel %d', ...
                continuous_config(u).electrode,continuous_config(u).channel);
          else
            continue;
          end
        end
        try
          reference = read_continuous(sprintf( ...
              '%s/%s/continuous/%s_day%d_electrode%02d(%03d)_channel%d.eeg', ...
              path_prefix,this_subject,this_subject,this_day, ...
              continuous_config(u).electrode,continuous_config(u).depth, ...
              continuous_config(u).channel),continuous_config(u), ...
              session([session(:).day] == this_day));
        catch
          error('read_continuous failed for reference channel of %s',filename);
        end
        % check that reference has same matching timestamps
        if ~all(arrayfun(@(s1,s2) numel(s1.timestamp) == numel(s2.timestamp), ...
            continuous,reference))
          error('reference has different number of samples from continuous');
        end
        if ~all(arrayfun(@(s1,s2) all(s1.timestamp == s2.timestamp), ...
            continuous,reference))
          error('reference timestamps do not match continuous timestamps');
        end
        % add reference so that continuous samples are now referenced to
        % ground
        for i = 1:numel(continuous)
          continuous(i).samples = continuous(i).samples + reference(i).samples;
          continuous(i).reference = [0 0];
        end
      end

      % append to end of lfp
      if (exist('lfp') == 1)
        lfp = [lfp; continuous];
      else
        lfp = continuous;
      end
    end

    % save continuous data for this subject/day
    try
      %{
      save(sprintf(['%s/%s/continuous/%s_day%d_lfp.mat'], ...
          path_prefix,this_subject,this_subject,this_day),'lfp');
      %}
    catch
      error('could not save lfp data for %s day %d ', ...
          this_subject,this_day);
    end
    %{
    try
      save(sprintf('%s/%s/continuous/%s_day%d_emg.mat', ...
          path_prefix,this_subject,this_subject, ...
          this_day),'emg');
    catch
      error('could not save EMG data for %s day %d', ...
          this_subject,this_day);
    end
    %}
  end
end

