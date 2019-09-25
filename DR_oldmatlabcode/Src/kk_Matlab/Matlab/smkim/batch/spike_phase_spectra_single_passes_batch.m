
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48', 'S58', 'S59', 'S60', 'S61'};
% List of the electrodes to use for LFP signal
lfp_electrodes = struct( ...
    'left'  , { 6,  4,  3,  5,  1}, ...
    'right' , {10,  9, 12,  8, 10} );

load('/data14/smkim/theta_state.mat');
load('/data14/smkim/linearized_passes.mat');

TS_PER_SEC = 1e4;

% We include a pass only if the overall mean running speed satisfies criterion
% and theta power ratio exceeds criterion for at least 90% of the pass
% duration. These are the same criteria as I use for receptive field estimation.
MINIMUM_OVERALL_RUNNING_SPEED = 10;
MINIMUM_THETA_STATE_FRACTION = 0.9;

% multitaper spectral estimation parameters
params = struct( ...
    'taper_family'    , 'slepian', ...
    'W'               , 0.05, ...
    'bins_per_cycle'  , 180, ...
    'nfft'            , 2^16, ...
    'frequency_band'  , [0.5 1.5], ...
    'num_jitters'     , 500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)

  files_list = dir(sprintf('%s/%s/spike/*_unit.mat',path_prefix,subjects{s}));
  for i = 1:numel(files_list)

    % Load single-unit spike data
    filename = sprintf('%s/%s/spike/%s', ...
        path_prefix,subjects{s},files_list(i).name);
    load(filename);
    disp(filename);

    % Remove epochs that aren't run epochs
    unit = unit(strncmp('run',{unit(:).epoch},3));

    % Sort by epoch
    [label, epochs] = grp2idx({unit(:).epoch});

    for j = 1:numel(epochs)
      % Grab indices of elements in unit that correspond to epochs{j}
      unit_idx = find(label == j);
      assert(all(strcmp(epochs{j},{unit(unit_idx).epoch})));

      for k = 1:numel(unit_idx)

        % Get matching theta-filtered LFP (recorded in same hemisphere)
        if (exist('theta') ~= 1) || ...
            ~any(struct_cmp(unit(unit_idx(k)),theta, ...
            {'subject','day','epoch','hemisphere'}))
          load(sprintf('%s/%s/continuous/%s_day%d_theta.mat', ...
              path_prefix,subjects{s},subjects{s},unit(unit_idx(k)).day));
        end
        theta_idx = find(struct_cmp(unit(unit_idx(k)),theta, ...
            {'subject','day','epoch','hemisphere'}));
        assert(isscalar(theta_idx));

        % Get matching theta state
        theta_state_idx = find(struct_cmp(unit(unit_idx(k)), ...
            theta_state,{'subject','day','epoch'}));
        assert(isscalar(theta_state_idx));
        theta_timerange = theta_state(theta_state_idx).timerange;

        % Get linearized passes that occur during this recording session and
        % satisfy overall mean speed criterion
        passes = linearized_passes( ...
            struct_cmp(unit(unit_idx(k)), ...
            linearized_passes,{'subject','day','epoch'}) & ...
            ([linearized_passes(:).overall_mean_speed]' >= ...
            MINIMUM_OVERALL_RUNNING_SPEED));

        % For each pass, we check the fraction of time during which theta power
        % exceeds threshold and include only if this theta fraction is above
        % criterion
        passes = passes( ...
            arrayfun(@(x) double(sum(length_intervals( ...
            intersect_intervals(x.timerange,theta_timerange)))) / ...
            double(sum(length_intervals(x.timerange))) >= ...
            MINIMUM_THETA_STATE_FRACTION,passes));

        % Collapse passes into a scalar struct and retain only metadata fields
        % of interest
        passes = struct( ...
            'timerange', {{passes(:).timerange}'}, ...
            'direction', {[passes(:).direction]'}, ...
            'overall_mean_speed', {[passes(:).overall_mean_speed]'} );

        % Compute spectrograms
        spike_phase_spectra(k,1) = windowed_spike_phase_spectra( ...
            unit(unit_idx(k)),theta(theta_idx),passes,params);
      
        % Clear passes for safety
        clear('passes');

      end
        
      % Save spectra for this subject/day/epoch
      outfilename = sprintf('%s/%s/%s', ...
          path_prefix,subjects{s}, ...
          [files_list(i).name(1:end-8) epochs{j} '_spike_phase_spectra.mat']);
      disp(outfilename);
      disp(spike_phase_spectra);
      save(outfilename,'spike_phase_spectra','-v7.3');

      clear('unit_idx','theta','spike_phase_spectra','outfilename');
    end
    clear('unit','label','epochs');
  end
end


