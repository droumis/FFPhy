

% also edit: is_continuous, is_spike, is_continuous_config, is_spike_config

% *_continuous.mat and *_emg.mat
% *_spike.mat
% *_unit.mat
% *_recording_sites.mat
% *_continuous_config.mat
% *_spike_config.mat


clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48','S58','S59','S60','S61'};
left_elects = 1:7;
right_elects = 8:14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)

  % recording_sites struct array
  load(sprintf('%s/%s/%s_recording_sites.mat', ...
      path_prefix,subjects{s},subjects{s}));
  [recording_sites(:).hemisphere] = deal('NA');
  [recording_sites( ...
      ismember([recording_sites(:).electrode],left_elects)).hemisphere] = ...
      deal('left');
  [recording_sites( ...
      ismember([recording_sites(:).electrode],right_elects)).hemisphere] = ...
      deal('right');
  assert(is_recording_sites(recording_sites));
  %
  save(sprintf('%s/%s/%s_recording_sites.mat', ...
      path_prefix,subjects{s},subjects{s}),'recording_sites');
  %

  % continuous_config struct array
  load(sprintf('%s/%s/continuous/%s_continuous_config.mat', ...
      path_prefix,subjects{s},subjects{s}));
  [continuous_config(:).hemisphere] = deal('NA');
  [continuous_config( ...
      ismember([continuous_config(:).electrode],left_elects)).hemisphere] = ...
      deal('left');
  [continuous_config( ...
      ismember([continuous_config(:).electrode],right_elects)).hemisphere] = ...
      deal('right');
  assert(is_continuous_config(continuous_config));
  %  
  save(sprintf('%s/%s/continuous/%s_continuous_config.mat', ...
      path_prefix,subjects{s},subjects{s}),'continuous_config');
  %

  % spike_config struct array
  load(sprintf('%s/%s/spike/%s_spike_config.mat', ...
      path_prefix,subjects{s},subjects{s}));
  [spike_config( ...
      ismember([spike_config(:).tetrode],left_elects)).hemisphere] = ...
      deal('left');
  [spike_config( ...
      ismember([spike_config(:).tetrode],right_elects)).hemisphere] = ...
      deal('right');
  assert(is_spike_config(spike_config));
  %
  save(sprintf('%s/%s/spike/%s_spike_config.mat', ...
      path_prefix,subjects{s},subjects{s}),'spike_config');
  %
 
  % continuous data 
  cont_files = dir(sprintf('%s/%s/continuous/*_continuous.mat', ...
      path_prefix,subjects{s}));
  for c = 1:numel(cont_files)
    load(sprintf('%s/%s/continuous/%s', ...
        path_prefix,subjects{s},cont_files(c).name));
    [continuous( ...
        ismember([continuous(:).electrode],left_elects)).hemisphere] = ...
        deal('left');
    [continuous( ...
        ismember([continuous(:).electrode],right_elects)).hemisphere] = ...
        deal('right');
    assert(is_continuous(continuous));
    %
    save(sprintf('%s/%s/continuous/%s', ...
        path_prefix,subjects{s},cont_files(c).name),'continuous');
    %
  end

  % emg data
  emg_files = dir(sprintf('%s/%s/continuous/*_emg.mat', ...
      path_prefix,subjects{s}));
  for e = 1:numel(emg_files)
    load(sprintf('%s/%s/continuous/%s', ...
        path_prefix,subjects{s},emg_files(e).name));
    [continuous(:).hemisphere] = deal('NA');
    assert(is_continuous(continuous));
    %  
    save(sprintf('%s/%s/continuous/%s', ...
        path_prefix,subjects{s},emg_files(e).name),'continuous');
    %
  end

  % spike data
  tet_folders = dir(sprintf('%s/%s/spike/tetrode*', ...
      path_prefix,subjects{s}));
  for t = 1:numel(tet_folders)
    spike_files = dir(sprintf('%s/%s/spike/%s/*_spike.mat', ...
      path_prefix,subjects{s},tet_folders(t).name));
    for w = 1:numel(spike_files)
      %{
      load(sprintf('%s/%s/spike/%s/%s', ...
          path_prefix,subjects{s},tet_folders(t).name,spike_files(w).name));
      %}
      spike = load(sprintf('%s/%s/spike/%s/%s', ...
          path_prefix,subjects{s},tet_folders(t).name,spike_files(w).name));
      [spike( ...
          ismember([spike(:).tetrode],left_elects)).hemisphere] = ...
          deal('left');
      [spike( ...
          ismember([spike(:).tetrode],right_elects)).hemisphere] = ...
          deal('right');
      assert(is_spike(spike));
      %
      save(sprintf('%s/%s/spike/%s/%s', ...
          path_prefix,subjects{s},tet_folders(t).name,spike_files(w).name), ...
          'spike');
      %
    end
  end

  % unit data
  unit_files = dir(sprintf('%s/%s/spike/*_unit.mat', ...
      path_prefix,subjects{s}));
  for u = 1:numel(unit_files)
    load(sprintf('%s/%s/spike/%s', ...
        path_prefix,subjects{s},unit_files(u).name));
    [unit( ...
        ismember([unit(:).tetrode],left_elects)).hemisphere] = ...
        deal('left');
    [unit( ...
        ismember([unit(:).tetrode],right_elects)).hemisphere] = ...
        deal('right');
    %
    save(sprintf('%s/%s/spike/%s', ...
        path_prefix,subjects{s},unit_files(u).name),'unit');
    %
  end

end

