%
load('/data14/smkim/normalized_power_iter.mat');
load('/data14/smkim/unit_iter.mat');

subjects = {'S48','S58','S59','S60','S61'};

uids = unique({unit_iter.uid});
%


TS_PER_SEC = 1e4;
POWER_RATIO_THRESHOLD = 2; % 2 deciBels
% Histogram bin edges for log10 ISI 
log10_isi_grid = (-3:0.05:1)';
bin_centers = 0.5*(log10_isi_grid(1:end-1) + log10_isi_grid(2:end));
clear('isi_summary');
isi_summary = struct( ...
    'uid'                   , {}, ...
    'region'                , {}, ...
    'theta_firing_rate'     , {}, ...
    'nontheta_firing_rate'  , {}, ...
    'theta_isi'             , {}, ...
    'nontheta_isi'          , {}, ...
    'theta_burst_index'     , {}, ...
    'nontheta_burst_index'  , {} );

for i = 1:numel(uids)
  unit = select(unit_iter,@(x) strcmp(x.uid,uids{i}));

  % total number of spikes
  theta_num_spikes = 0;
  % number of single-spike (non-burst) events
  theta_num_single_spikes = 0; 
  % number of burst events; this is NOT the same as the number of spikes that
  % occurred during bursts
  theta_num_bursts = 0; 
  % amount of time spend in theta state
  theta_intervals_duration = 0; % in seconds
  % vector of ISI bin counts
  theta_isi_counts = zeros(size(log10_isi_grid));

  % total number of spikes
  nontheta_num_spikes = 0;
  % number of single-spike (non-burst) events
  nontheta_num_single_spikes = 0; 
  % number of burst events; this is NOT the same as the number of spikes that
  % occurred during bursts
  nontheta_num_bursts = 0; 
  % amount of time spend in theta state
  nontheta_intervals_duration = 0; % in seconds
  % vector of ISI bin counts
  nontheta_isi_counts = zeros(size(log10_isi_grid));

  for j = 1:numel(unit)
    normalized_power = select(normalized_power_iter,@(x) struct_cmp( ...
        unit(j),x,{'subject','day','epoch'}));
    % Overall time interval (in timestamp units)
    overall_interval = double(unit(j).timerange);
    % Find theta time intervals (in timestamp units)
    theta_intervals = find_intervals(double(normalized_power.timestamp), ...
        normalized_power.theta_spindle_power_ratio, ...
        @(x) x > POWER_RATIO_THRESHOLD,1);
    % Non-theta intervals are the complement (in timestamp units)
    nontheta_intervals = diff_intervals(overall_interval,theta_intervals);
    % Compute lengths of the intervals
    theta_intervals_duration = theta_intervals_duration + ...
        length_intervals(theta_intervals)/TS_PER_SEC;
    nontheta_intervals_duration = nontheta_intervals_duration + ...
        length_intervals(nontheta_intervals)/TS_PER_SEC;
    % Select spikes that occur in theta versus nontheta intervals
    theta_spikes = ismember_intervals(double(unit(j).timestamp), ...
        theta_intervals);
    nontheta_spikes = ismember_intervals(double(unit(j).timestamp), ...
        nontheta_intervals);
    theta_num_spikes = theta_num_spikes + nnz(theta_spikes);
    nontheta_num_spikes = nontheta_num_spikes + nnz(nontheta_spikes);
    % Find single-spike events and burst events
    if any(theta_spikes)
      theta_num_single_spikes = theta_num_single_spikes + ...
          nnz(unit(j).burst(2).flag(theta_spikes) == 0);
      theta_num_bursts = theta_num_bursts + ...
          nnz(unit(j).burst(2).flag(theta_spikes) == 1);
    end
    if any(nontheta_spikes)
      nontheta_num_single_spikes = nontheta_num_single_spikes + ...
          nnz(unit(j).burst(2).flag(nontheta_spikes) == 0);
      nontheta_num_bursts = nontheta_num_bursts + ...
          nnz(unit(j).burst(2).flag(nontheta_spikes) == 1);
    end
    % Histogram counts of log10 interspike intervals
    if any(theta_spikes)
      tmp = histc(log10(unit(j).isi_post(theta_spikes)),log10_isi_grid);
      if (size(tmp,1) < size(tmp,2))
        tmp = tmp';
      end
      theta_isi_counts = theta_isi_counts + tmp;
    end
    if any(nontheta_spikes)
      tmp = histc(log10(unit(j).isi_post(nontheta_spikes)),log10_isi_grid);
      if (size(tmp,1) < size(tmp,2))
        tmp = tmp';
      end
      nontheta_isi_counts = nontheta_isi_counts + tmp;
    end
  end  
  % Append to the isi_summary
  isi_summary(end+1,1) = struct( ...
      'uid'                   , unit(j).uid           , ...
      'region'                , unit(j).region        , ...
      'theta_firing_rate'     , theta_num_spikes / theta_intervals_duration , ...
      'nontheta_firing_rate'  , nontheta_num_spikes / nontheta_intervals_duration , ...
      'theta_isi'             , theta_isi_counts      , ...
      'nontheta_isi'          , nontheta_isi_counts   , ...
      'theta_burst_index'     , theta_num_bursts/(theta_num_bursts+theta_num_single_spikes), ...
      'nontheta_burst_index'  , nontheta_num_bursts/(nontheta_num_bursts+nontheta_num_single_spikes) );
end



