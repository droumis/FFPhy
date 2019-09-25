
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

% Parameters for spacing windows, in clock time (seconds, not theta cycles). We
% want uniform dense coverage over the times of interest, but if we make the
% window spacing too fine then the computation takes too long.
step_size = 3;
% Width of window in seconds.
window_size = 6;

% multitaper spectral estimation parameters
params = struct( ...
    'taper_family'    , 'slepian', ...
    'W'               , 0.05, ...
    'bins_per_cycle'  , 360, ...
    'nfft'            , 2^16, ...
    'frequency_band'  , [0 2], ...
    'num_jitters'     , 200);

% Slepian tapers are optimized for spectral concentration (sidelobe
% suppression), whereas sinusoidal tapers minimize smoothing bias (Gibbs
% ringing) within local bandwidth.

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
        % Load matching theta and select desired LFP reference electrode
        if (exist('theta') ~= 1) || ...
            ~any(struct_cmp(unit(unit_idx(k)),theta, ...
            {'subject','day','epoch','environment','hemisphere'}))
          load(sprintf('%s/%s/continuous/%s_day%d_theta.mat', ...
              path_prefix,subjects{s},subjects{s},unit(unit_idx(k)).day));
        end
        theta_idx = find(struct_cmp(unit(unit_idx(k)),theta, ...
            {'subject','day','epoch','environment','hemisphere'}));

        % TODO:
        % Evenly spaced overlapping windows to span the data timerange, with
        % margins at the start and end of the session
        duration = double(unit(unit_idx(k)).timerange(end) - ...
            unit(unit_idx(k)).timerange(1))/TS_PER_SEC;
        num_windows = floor(1 + (duration-window_size)/step_size);
        windows.timerange = mat2cell(bsxfun(@plus, ...
          unit(unit_idx(k)).timerange(1) + ...
          uint32(TS_PER_SEC*(rem(duration-window_size,step_size)/2 + ...
          step_size*(0:(num_windows-1))')), ...
          uint32([0 window_size]*TS_PER_SEC)),ones([num_windows 1]),[2]);

        %{
        % Include only windows that satisfy some state criterion
        timerange = intersect_intervals( ...
            cell2mat(linearized_passes( ...
            struct_cmp(unit(unit_idx(k)),linearized_passes, ...
            {'subject','day','epoch','environment'})).timerange), ...
            theta_state(struct_cmp(unit(unit_idx(k)),theta_state, ...
            {'subject','day','epoch','environment'})).timerange);
        %}

        % Compute spectrograms
        spike_phase_spectra(k,1) = windowed_spike_phase_spectra( ...
            unit(unit_idx(k)),theta(theta_idx),windows,params);
      
        % Clear windows struct for safety
        clear('windows');
      end
        
      % DEBUG
      return;

      % Save spectra for this subject/day/epoch
      outfilename = sprintf('%s/%s/%s', ...
          path_prefix,subjects{s}, ...
          [files_list(i).name(1:end-8) epochs{j} '_spike_phase_spectra.mat']);
      disp(outfilename);
      disp(spike_phase_spectra);
      save(outfilename,'spike_phase_spectra','-v7.3');

      clear('unit_idx','spike_phase_spectra','outfilename');
    end
    clear('unit','label','epochs');
  end
end

