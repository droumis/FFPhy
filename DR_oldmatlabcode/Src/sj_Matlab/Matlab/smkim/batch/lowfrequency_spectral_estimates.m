
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
tgrid = (1:samples_per_window)'/Fs - mean(1:samples_per_window)/Fs;
% Construct tapers on grid
% For a given T*W product, we get spectral concentration within the frequency
% interval [-W, +W]
TW = 2; 
[tapers, eigs] = dpsschk([TW 2*TW-1],samples_per_window,Fs);
% Number of points in fft of prolates
nfft = 2^14;
% Indices to select desired components out of the entire frequency range
fpass = [1, 30];
[fgrid, fidx] = getfgrid(Fs,nfft,fpass);
% window_overlap is specified as a fraction of diff(window_size)
window_overlap = 0.8;
stepsize = diff(window_size) * (1 - window_overlap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
% estimate low-frequency spectra for LFP
% Load each *_lfp.mat file in turn and process all electrodes
for s = 1:numel(subjects)

  % Load prewhitening filters for this subject
  load(sprintf('%s/%s/continuous/%s_lfp_prewhitening_filters.mat', ...
      path_prefix,subjects{s},subjects{s}));

  files_list = dir(sprintf('%s/%s/continuous/*_lfp.mat',path_prefix,subjects{s}));
  for i = 1:numel(files_list)

    lfp_spectra = [];

    filename = sprintf('%s/%s/continuous/%s', ...
        path_prefix,subjects{s},files_list(i).name);
    load(filename);

    % Select desired LFP reference electrodes
    lfp = lfp(([lfp(:).electrode] == lfp_electrodes(s).left) | ...
        ([lfp(:).electrode] == lfp_electrodes(s).right));

    for j = 1:numel(lfp)
      % Whiten using the prewhitening filter
      Hd = prewhitening_filters(struct_cmp(lfp(j), ...
          prewhitening_filters,{'subject','day','electrode'})).filter_object;
      % whiten_continuous uses filtfilthd
      whitened = filter_continuous(lfp(j),Hd);
      assert(isscalar(whitened));
      assert(isequal(whitened.timestamp,lfp(j).timestamp));
      
      lfp_spectra(j,1).subject = whitened.subject;
      lfp_spectra(j).day = whitened.day;
      lfp_spectra(j).epoch = whitened.epoch;
      lfp_spectra(j).environment = whitened.environment;
      lfp_spectra(j).region = whitened.region;
      lfp_spectra(j).hemisphere = whitened.hemisphere;
      lfp_spectra(j).depth = whitened.depth;
      lfp_spectra(j).electrode = whitened.electrode;
      lfp_spectra(j).channel = whitened.channel;
      lfp_spectra(j).reference = whitened.reference;
      lfp_spectra(j).Fs = whitened.Fs;
      % Append name of the *_lfp.mat file to the sources cell array
      lfp_spectra(j).sources = [filename, whitened.sources];
      % Add prewhitening filter object
      lfp_spectra(j).prewhitening_filter = Hd;

      assert(numel(lfp_spectra) == j);
      
      % Divide the duration of LFP record into overlapping windows,
      % discarding fractional-window leftovers at the start and end
      n_windows = 1 + floor((double(whitened.timestamp(end) - ...
          whitened.timestamp(1))/TS_PER_SEC - diff(window_size)) / stepsize);
      assert(n_windows > 0);
      bincenters = double(whitened.timestamp(1))/TS_PER_SEC + ...
          rem(double(whitened.timestamp(end) - ...
          whitened.timestamp(1))/TS_PER_SEC,stepsize)/2 + ...
          diff(window_size)/2 - stepsize + stepsize*(1:n_windows)';

      % Grab LFP samples in each window
      data = zeros([samples_per_window, n_windows]);
      for k = 1:n_windows
        samples_idx = [ ...
            find(double(whitened.timestamp)/TS_PER_SEC < bincenters(k), ...
            floor(samples_per_window/2),'last');
            find(double(whitened.timestamp)/TS_PER_SEC >= bincenters(k), ...
            ceil(samples_per_window/2),'first') ];
          assert(numel(samples_idx) == samples_per_window);
          data(:,k) = double(whitened.samples(samples_idx));
      end

      % Compute sliding-window estimates of spectral density for whitened signal
      lfp_spectra(j).timestamp = uint32(TS_PER_SEC*bincenters);
      J_w = mtfftc(data,tapers,nfft,Fs);
      J_w = J_w(fidx,:,:);
      % Here we simply average the projections instead of using Thomson's
      % adaptive method. As long as we are using the (2*TW - 1) tapers with
      % eigenvalues close to unity, the average is very close to optimal
      lfp_spectra(j).normalized_spectral_density = ...
          squeeze(mean(conj(J_w) .* J_w,2))';

      % confidence interval
      %Serr=specerr(S,J,err,trialave,Nsp);
      % Inputs:
      % S - spectrum
      % J - tapered fourier transforms 
      % err - [errtype p] (errtype=1 - asymptotic estimates; errchk=2 - Jackknife estimates; 
      %                   p - p value for error estimates)
      % trialave - 0: no averaging over trials/channels
      %            1 : perform trial averaging
      % numsp    - number of spikes in each channel. specify only when finite
      %            size correction required (and of course, only for point
      %            process data)
      %
      clear('J_w','data','whitened');      

      % Compute sliding-window estimates of spectral density for raw signal
      data = zeros([samples_per_window, n_windows]);
      for k = 1:n_windows
        samples_idx = [ ...
            find(double(lfp(j).timestamp)/TS_PER_SEC < bincenters(k), ...
            floor(samples_per_window/2),'last');
            find(double(lfp(j).timestamp)/TS_PER_SEC >= bincenters(k), ...
            ceil(samples_per_window/2),'first') ];
          assert(numel(samples_idx) == samples_per_window);
          data(:,k) = double(lfp(j).samples(samples_idx));
      end
      J_l = mtfftc(data,tapers,nfft,Fs);
      J_l = J_l(fidx,:,:);
      lfp_spectra(j).spectral_density = squeeze(mean(conj(J_l) .* J_l,2))';
      clear('J_l','data');

    end

    % Deal the common meta-data
    [lfp_spectra(:).frequency] = deal(fgrid);
    [lfp_spectra(:).window_size] = deal(window_size);
    [lfp_spectra(:).window_overlap] = deal(window_overlap);
    [lfp_spectra(:).TW] = deal(TW);
    [lfp_spectra(:).num_tapers] = deal(size(tapers,2));

    % Save
    disp(lfp_spectra);
    outfilename = [filename(1:end-4) '_lowfrequencyspectra.mat'];
    disp(outfilename);
    save(outfilename,'lfp_spectra');
    clear('lfp','lfp_spectra','outfilename');
  end

end
%}

%{
% spectral estimates for single-unit spike data
% Load each *_unit.mat file in turn and process all electrodes
for s = 1:numel(subjects)
  files_list = dir(sprintf('%s/%s/spike/*_unit.mat',path_prefix,subjects{s}));
  for i = 1:numel(files_list)
    filename = sprintf('%s/%s/spike/%s', ...
        path_prefix,subjects{s},files_list(i).name);
    disp(filename);
    load(filename);

    % split up each unit data struct into epochs
    unique_epochs = unique({unit(:).epoch});
    for e = 1:numel(unique_epochs)
      epoch = unique_epochs{e};
      epoch_lookup = find(strcmp(epoch,{unit(:).epoch}));
      unit_spectra = rmfield(unit(epoch_lookup), ...
          {'passbands','thresholds','Fs', ...
          'samples_before_trigger','amplitude_cutoff','samples', ...
          'maximum_amplitude_channel','shape','isi_pre','isi_post','burst'});

      for j = 1:numel(unit_spectra)
        % Divide the duration of spike record into overlapping windows,
        % discarding fractional-window leftovers at the start and end
        n_windows = 1 + floor((double(unit(epoch_lookup(j)).timerange(end) - ...
            unit(epoch_lookup(j)).timerange(1))/TS_PER_SEC - ...
            diff(window_size)) / stepsize);
        assert(n_windows > 0);
        bincenters = double(unit(epoch_lookup(j)).timerange(1))/TS_PER_SEC + ...
            rem(double(unit(epoch_lookup(j)).timerange(end) - ...
            unit(epoch_lookup(j)).timerange(1))/TS_PER_SEC,stepsize)/2 + ...
            diff(window_size)/2 - stepsize + stepsize*(1:n_windows)';

        % Grab spikes in each window
        data = struct('times',cell([n_windows, 1]));
        spiketimes = double(unit(epoch_lookup(j)).timestamp)/TS_PER_SEC;
        for k = 1:n_windows
          unit_idx = find( (spiketimes - bincenters(k) >= min(tgrid)) & ...
              (spiketimes - bincenters(k) <= max(tgrid)) );
          % Collect spike times relative to the bin center (this makes it easier
          % to vectorize the projections onto the tapers)
          data(k).times = spiketimes(unit_idx) - bincenters(k);
        end

        % Compute spectrum in each window
        unit_spectra(j).timestamp = uint32(TS_PER_SEC*bincenters);
        [J_unit, junk, num_spikes] = mtfftpt(data,tapers,nfft,tgrid,fgrid,fidx);
        unit_spectra(j).spectral_density = squeeze(mean(conj(J_unit) .* J_unit,2))';
        % Number of spikes that contributed to the estimate in each window
        unit_spectra(j).num_spikes = num_spikes(:);

        % confidence interval
        %Serr=specerr(S,J,err,trialave,Nsp);
        % Inputs:
        % S - spectrum
        % J - tapered fourier transforms 
        % err - [errtype p] (errtype=1 - asymptotic estimates; errchk=2 - Jackknife estimates; 
        %                   p - p value for error estimates)
        % trialave - 0: no averaging over trials/channels
        %            1 : perform trial averaging
        % numsp    - number of spikes in each channel. specify only when finite
        %            size correction required (and of course, only for point
        %            process data)
        %

        % Append name of the *_unit.mat file to the sources cell array
        unit_spectra(j).sources = [filename, unit(epoch_lookup(j)).sources];
      end
    
      % Deal the common meta-data
      [unit_spectra(:).window_size] = deal(window_size);
      [unit_spectra(:).window_overlap] = deal(window_overlap);
      [unit_spectra(:).frequency] = deal(fgrid);
      [unit_spectra(:).TW] = deal(TW);
      [unit_spectra(:).num_tapers] = deal(size(tapers,2));
      % Fs is sampling resolution of the tapers, not of the spikes
      [unit_spectra(:).Fs] = deal(Fs);

      % Save
      disp(unit_spectra);
      outfilename = [filename(1:end-8) epoch '_unit_spectra.mat'];
      disp(outfilename);
      save(outfilename,'unit_spectra');
      clear('unit_spectra','epoch','outfilename');
    end
    clear('unit');
  end
end
%}

%

