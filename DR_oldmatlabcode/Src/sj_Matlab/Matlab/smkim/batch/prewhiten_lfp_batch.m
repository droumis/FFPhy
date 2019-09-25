
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48', 'S58', 'S59', 'S60', 'S61'};
% List of the electrodes to use for LFP signal
lfp_electrodes = struct( ...
    'left'  , { 6,  4,  3,  5,  1}, ...
    'right' , {10,  9, 12,  8, 10} );

TS_PER_SEC = 1e4;
% LFP is sampled at 1500 samples/second; we bin the spike times at the same
% resolution
Fs = 1500;
% Size of the sliding window
window_size = [-1, +1];
samples_per_window = round(Fs*diff(window_size));
% Timegrid at which the tapers are sampled, measured relative to window center
% in the same units as 1/Fs
t = (1:samples_per_window)'/Fs - mean(1:samples_per_window)/Fs;
% Construct tapers on grid
% For a given T*W product, we get spectral concentration within the frequency
% interval [-W, +W]
TW = 2; 
[tapers, eigs] = dpsschk([TW 2*TW-1],samples_per_window,Fs);
% Number of points in fft of prolates
nfft = 2^14;
% Indices to select desired components out of the entire frequency range; here
% we use the entire passband of the NSpike analog filters
fpass = [1, 400];

[fgrid, fidx] = getfgrid(Fs,nfft,fpass);
% window_overlap is specified as a fraction of diff(window_size)
window_overlap = 0.5;
stepsize = diff(window_size) * (1 - window_overlap);

% regularization penalty for smoothing spline
smoothing_parameter = 0.1;
% to avoid theta and its first harmonic, we exclude 4-18 Hz from fitting. The
% smoothing spline guarantees that the fit will be well-behaved over this
% interval
fit_frequencies = (fgrid < 4) | (fgrid > 18);

% filter length
n_filter = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:numel(subjects)

  prewhitening_filters = struct( ...
      'subject'         , {}, ...
      'day'             , {}, ...
      'electrode'       , {}, ...
      'hemisphere'      , {}, ...
      'log10frequency_dBpower_spline_fit', {}, ...
      'filter_object'   , {} );

  % Load all position data
  load(sprintf('%s/%s/behavior/%s_position.mat', ...
      path_prefix,subjects{i},subjects{i}));

  % For each reference LFP electrode, collect all LFP data across epochs on each
  % day
  for day = 0:7

    include_intervals = {};
    S.left = [];
    S.right = [];
    % use only the rest sessions, because LFP recorded during run sessions is
    % corrupted by low-frequency movement artifact.
    files_list = dir(sprintf('/%s/%s/continuous/%s_day%d_rest*_lfp.mat', ...
        path_prefix,subjects{i},subjects{i},day));
    if isempty(files_list)
      continue;
    end

    for k = 1:numel(files_list)
      load(sprintf('%s/%s/continuous/%s', ...
          path_prefix,subjects{i},files_list(k).name));

      % Find matching time intervals during which rat is not immobile
      pos_idx = find(struct_cmp(lfp(1),position,{'subject','day','epoch'}));
      include_intervals = find_intervals( ...
          position(pos_idx).timestamp,double(position(pos_idx).stopped), ...
          @(x) x == 0,uint32(10));

      % Look up the designated electrode in each hemisphere
      lfp_idx.left = find( ...
            ([lfp(:).electrode] == lfp_electrodes(i).left) & ...
            strcmp('box',{lfp(:).environment}));
      lfp_idx.right = find( ...
            ([lfp(:).electrode] == lfp_electrodes(i).right) & ...
            strcmp('box',{lfp(:).environment}));
      assert(isscalar(lfp_idx.left) && isscalar(lfp_idx.right)); 

      % Determine time bins
      common_timerange = [ ...
          max(lfp(lfp_idx.left).timestamp(1), ...
          lfp(lfp_idx.right).timestamp(1)), ...
          min(lfp(lfp_idx.left).timestamp(end), ...
          lfp(lfp_idx.right).timestamp(end)) ];

      % Divide the duration of LFP record into overlapping windows,
      % discarding fractional-window leftovers at the start and end
      n_windows = 1 + floor((double(common_timerange(end) - ...
          common_timerange(1))/TS_PER_SEC - diff(window_size)) / stepsize);
      assert(n_windows > 0);
      bincenters = double(common_timerange(1))/TS_PER_SEC + ...
          rem(double(common_timerange(end) - ...
          common_timerange(1))/TS_PER_SEC,stepsize)/2 + ...
          diff(window_size)/2 - stepsize + stepsize*(1:n_windows)';

      % Grab LFP samples in each window
      data.left = zeros([samples_per_window, n_windows]);
      data.right = zeros([samples_per_window, n_windows]);
      for k = 1:n_windows
        samples_idx.left = [ ...
            find(double(lfp(lfp_idx.left).timestamp)/TS_PER_SEC < ...
            bincenters(k),floor(samples_per_window/2),'last');
            find(double(lfp(lfp_idx.left).timestamp)/TS_PER_SEC >= ...
            bincenters(k),ceil(samples_per_window/2),'first') ];
        samples_idx.right = [ ...
            find(double(lfp(lfp_idx.right).timestamp)/TS_PER_SEC < ...
            bincenters(k),floor(samples_per_window/2),'last');
            find(double(lfp(lfp_idx.right).timestamp)/TS_PER_SEC >= ...
            bincenters(k),ceil(samples_per_window/2),'first') ];
        assert(isequal(numel(samples_idx.left),numel(samples_idx.right), ...
            samples_per_window));
        data.left(:,k) = ...
            double(lfp(lfp_idx.left).samples(samples_idx.left));
        data.right(:,k) = ...
            double(lfp(lfp_idx.right).samples(samples_idx.right));
      end

      % Compute sliding-window estimates of spectral density
      J_w.left = mtfftc(data.left,tapers,nfft,Fs);
      J_w.left = J_w.left(fidx,:,:);
      J_w.right = mtfftc(data.right,tapers,nfft,Fs);
      J_w.right = J_w.right(fidx,:,:);

      % Select those windows that fall within the desired intervals 
      include_idx = ismember_intervals(uint32(TS_PER_SEC*bincenters), ...
          include_intervals);
      J_w.left = J_w.left(:,:,include_idx);
      J_w.right = J_w.right(:,:,include_idx);

      % Here we simply average the projections instead of using Thomson's
      % adaptive method. As long as we are using the (2*TW - 1) tapers with
      % eigenvalues close to unity, the average is very close to optimal
      S.left = [ S.left; squeeze(mean(conj(J_w.left) .* J_w.left,2))' ];
      S.right = [ S.right; squeeze(mean(conj(J_w.right) .* J_w.right,2))' ];

      clear('J_w');
    end

    % Fit the median of log-log spectral density profile. We want frequency to also
    % be on a logarithmic scale because the fit needs to be more wiggly to fit the
    % low frequencies
    sp.left = spaps(log10(fgrid(fit_frequencies)), ...
        median(db(S.left(:,fit_frequencies),'power'),1),smoothing_parameter);
    sp.right = spaps(log10(fgrid(fit_frequencies)), ...
        median(db(S.right(:,fit_frequencies),'power'),1),smoothing_parameter);

    % we want the filter gain to be inversely related to the spectral power
    % profile, and we want it to be scaled so that when the filter is applied
    % once forwards and once backwards it yields the desired "white" profile
    h.left = -0.25*fnval(sp.left,log10(fgrid));
    h.right = -0.25*fnval(sp.right,log10(fgrid));
    % normalize so that power is approximately preserved in theta band; this is
    % mainly for aesthetic reasons but also to keep values in a reasonable
    % floating-point range
    h.left = h.left - mean(h.left((fgrid > 5) & (fgrid < 10)));
    h.right = h.right - mean(h.right((fgrid > 5) & (fgrid < 10)));
    % Extrapolate to avoid edge effects
    h.left = [2*h.left(1)-h.left(2), h.left, 2*h.left(end)-h.left(end-1)];
    h.right = [2*h.right(1)-h.right(2), h.right, 2*h.right(end)-h.right(end-1)];
    % convert from deciBels to linear gains
    h.left = 10.^(h.left/10);
    h.right = 10.^(h.right/10);
    % construct prewhitening filter
    d.left = fdesign.arbmag(n_filter,[0, 2*fgrid/Fs, 1],h.left);
    d.right = fdesign.arbmag(n_filter,[0, 2*fgrid/Fs, 1],h.right);
    Hd.left = design(d.left,'fir');
    Hd.right = design(d.right,'fir');

    prewhitening_filters = [ prewhitening_filters; struct( ...
        'subject'         , {subjects{i}; subjects{i}}, ...
        'day'             , {day; day}, ...
        'electrode'       , {lfp_electrodes(i).left; lfp_electrodes(i).right}, ...
        'hemisphere'      , {'left'; 'right'}, ...
        'log10frequency_dBpower_spline_fit', {sp.left; sp.right}, ...
        'filter_object'   , {Hd.left; Hd.right} ) ];

    clear('S','sp','h','d','Hd');
  end
  disp(prewhitening_filters);
  save(sprintf('%s/%s/continuous/%s_lfp_prewhitening_filters.mat', ...
      path_prefix,subjects{i},subjects{i}),'prewhitening_filters');
end



