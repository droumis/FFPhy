function measures = receptive_field_measures(receptive_fields)
%RECEPTIVE_FIELD_MEASURES Various measures that are commonly used to quantify place fields.
%
%   MEASURES = RECEPTIVE_FIELD_MEASURES(RECEPTIVE_FIELDS) computes various
%   measures on the position-phase receptive fields in RECEPTIVE_FIELDS.
%   RECEPTIVE_FIELDS must be a 
%
%   References:
%   [1] Skaggs W.E., McNaughton B.L., Wilson M.A., Barnes C.A. (1996) Theta
%       phase precession in hippocampal neuronal populations and the compression
%       of temporal sequences. _Hippocampus_ 6:149-172.
%   [2] Battaglia F.P., Sutherland G.R., McNaughton B.L. (2004) Local sensory
%       cues and place cell directionality: additional evidence of prospective
%       coding in the hippocampus. _Journal of Neuroscience_ 24:4541-4550.
%   [3] Brun V.H., Solstad T., Kjelstrup K.B., Fyhn M., Witter M.P., Moser E.I.,
%       Moser M.B. (2008) Progressive increase in grid scale from dorsal to
%       ventral medial entorhinal cortex. _Hippocampus_ 18:1200-1212.
%
%Depends on:
%   NANSUM (MATLAB Statistics toolbox)
%   STRUCT_CMP (written by SMK)
%   IS_TIMERANGE (written by smk)
%   DIFF_INTERVALS (written by smk)
%   ISMEMBER_INTERVALS (written by smk)
%   DENSITY_MAP (written by smk)
%
%Written by SMK, 2009 November 10.
%

% activity fraction (sparsity)
% directionality
% field size (above some threshold)

  TS_PER_SEC = 1e4;

  if (exist('LINEARIZED_PLACE_FIELD') ~= 2)
    error(['LINEARIZED_PLAE_FIELD_MEASURES depends on m-file ' ...
        'LINEARIZED_PLACE_FIELD (written by SMK)']);
  end
  if (exist('NANSUM') ~= 2)
    error(['LINEARIZED_PLAE_FIELD_MEASURES depends on m-file ' ...
        'NANSUM (MATLAB Statistics Toolbox)']);
  end

  % Check that rightbound_timerange and leftbound_timerange don't overlap
  if is_timerange(rightbound_timerange) 
    rightbound_timerange = {rightbound_timerange};
  end
  if is_timerange(leftbound_timerange) 
    leftbound_timerange = {leftbound_timerange};
  end
  if ~iscell(rightbound_timerange) || ~iscell(leftbound_timerange) || ...
      ~all(cellfun(@is_timerange,rightbound_timerange)) || ...
      ~all(cellfun(@is_timerange,leftbound_timerange)) || ...
      ~isequal(size(rightbound_timerange),size(leftbound_timerange))
      any(cellfun(@(c1,c2) ~isempty(intersect_intervals(c1,c2)), ...
      rightbound_timerange,leftbound_timerange))
    error(['RIGHTBOUND_TIMERANGES and LEFTBOUND_TIMERANGES must be ' ...
      'equal-sized cell arrays containing uint32 timestamps intervals, ' ...
      'and corresponding sets of intervals must not overlap']);
  end

  % I'm lazy: pass along everything to LINEARIZED_PLACE_FIELD to catch errors
  try
    measures.ratemap.rightbound = linearized_place_field(unit,linearized, ...
        rightbound_timerange,params);
    measures.ratemap.leftbound = linearized_place_field(unit,linearized, ...
        leftbound_timerange,params);
  catch
    error('Could not estimate linearized place field. Bad inputs?');
  end
  
  % Combine data from leftbound and rightbound
  firing_rate_per_bin = [ measures.ratemap.rightbound.rate(:); ...
      measures.ratemap.leftbound.rate(:) ];
  measures.spatially_averaged_mean_rate = mean(firing_rate_per_bin);
  measures.peak_rate = max(firing_rate_per_bin);

  % Compute spatial information per spike and sparsity according to (Skaggs et
  % al., 1996)
  total_spike_count = 0;
  total_duration = 0;
  for i = 1:numel(unit)
    combined_intervals = union_intervals( ...
        rightbound_timerange{i},leftbound_timerange{i});
    total_spike_count = total_spike_count + ...
        nnz(ismember_intervals(unit(i).timestamp,combined_intervals));
    total_duration = total_duration + ...
        sum(double(length_intervals(combined_intervals)))/TS_PER_SEC;
  end
  measures.overall_mean_rate = total_spike_count / total_duration;
  occupancy_probability = [ measures.ratemap.rightbound.occupancy(:); ...
      measures.ratemap.leftbound.occupancy(:) ];
  occupancy_probability = occupancy_probability / sum(occupancy_probability);
  measures.spatial_information_per_spike = nansum( ...
      occupancy_probability .* ...
      (firing_rate_per_bin/measures.overall_mean_rate) .* ...
      log2(firing_rate_per_bin/measures.overall_mean_rate) );
  measures.sparsity = ...
      sum(occupancy_probability .* firing_rate_per_bin)^2 / ...
      sum(occupancy_probability .* firing_rate_per_bin.^2);
  
  % Compute (bi)directionality measure (Battaglia et al., 2004)
  profile.rightbound = measures.ratemap.rightbound.rate(:) / ...
      sum(measures.ratemap.rightbound.rate);
  profile.leftbound = measures.ratemap.leftbound.rate(:) / ...
      sum(measures.ratemap.leftbound.rate);
  numbins = (numel(params.linearized_position_bin_edges) - 1);
  binsize = mean(diff(params.linearized_position_bin_edges));
  offsets = ( -round(numbins/2) : +round(numbins/2) )';
  assert(isvector(offsets));
  measures.bidirectionality.misalignment = binsize * offsets;
  measures.bidirectionality.overlap = zeros(size(offsets));
  for i = 1:numel(offsets)
    r = profile.rightbound((1 + offsets(i)*(offsets(i) > 0)): ...
        (numbins + offsets(i)*(offsets(i) < 0)));
    l = profile.leftbound((1 - offsets(i)*(offsets(i) < 0)): ...
        (numbins - offsets(i)*(offsets(i) > 0)));
    measures.bidirectionality.overlap(i) = 2*sum(min([r,l],[],2)) / sum(r + l);
  end

  % Compute spatial autocorrelation of firing rate along linear track (Brun et
  % al., 2004), combining rightbound and leftbound
  [c.rightbound, lags] = xcorr(measures.ratemap.rightbound.rate - ...
      measures.spatially_averaged_mean_rate,'biased');
  [c.leftbound, lags] = xcorr(measures.ratemap.leftbound.rate - ...
      measures.spatially_averaged_mean_rate,'biased');
  measures.periodicity.lags = binsize * lags;
  measures.periodicity.spatial_autocovariance = ...
      0.5*(c.rightbound + c.leftbound);

end % end main function



