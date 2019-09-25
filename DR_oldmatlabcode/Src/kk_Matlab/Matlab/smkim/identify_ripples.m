function ripplestruct = identify_ripples(ripple,threshold,min_separation,min_duration)
%IDENTIFY_RIPPLES Find ripple times by thresholding bandpass filtered data.
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This is a draft. Do not use for mission-critical applications!
%
%TODO: obtain peak, mean frequency during ripple
%
%
%   RIPPLESTRUCT = IDENTIFY_RIPPLES(RIPPLE,THRESHOLD,MIN_SEPARATION,MIN_DURATION)
%
%   fields: 
%     timerange (Nx2 uint32 matrix)
%     mean power (Nx1 vector)
%     mean frequency (Nx1 vector)
%
%   References:
%   
%   Lubenov & Siapas (2009) Hippocampal theta oscillations are travelling waves.
%   _Nature_ 459:534-539.
%
%   Wierzynski et al. (2009) State-dependent spike-timing relationships between
%   hippocampal and prefrontal circuits during sleep. _Neuron_ 61:587-596.
%
%Depends on:
%   FIND_INTERVALS (written by SMK)
%
%Written by SMK, 2009 November 18.
%

if (exist('find_intervals') ~= 2)
  error(['IDENTIFY_RIPPLES depends on the m-file FIND_INTERVALS ' ...
      '(written by SMK)']);
end

if ~is_continuous(ripple) || any(cellfun(@(c) ~isfloat(c) || isreal(c), ...
      {ripple.samples}))
  error('RIPPLE must be a continuous data struct with complex float samples field');
end

if ~isscalar(threshold) || ~isfloat(threshold) || ~isreal(threshold) || ...
    ~(threshold > 0)
  error('THRESHOLD must be a positive real floating-point scalar');
end

if ~isscalar(min_separation) || ~isreal(min_separation) || ...
    ~isa(min_separation,'uint32')
  error('MIN_SEPARATION must be a real uint32 scalar');
end

if ~isscalar(min_duration) || ~isreal(min_duration) || ...
    ~isa(min_duration,'uint32') || ~(min_duration > min_separation)
  error('MIN_DURATION must be a real uint32 scalar, greater than MIN_SEPARATION');
end

for i = 1:numel(ripple)
  
  % Find threshold crossings
  threshold_crossings = uint32(find_intervals( ...
      double(ripple(i).timestamp),abs(ripple(i).samples),@(x) x > threshold));
  
  % Combine threshold-crossings that are less than min_separation apart from
  % each other
  combine_flag = ...
      threshold_crossings(2:end,1) - threshold_crossings(1:end-1,2) < ...
      min_separation;
  candidate_events = zeros( ...
      [size(threshold_crossings,1)-nnz(combine_flag), 2],'uint32');
  % Each continuous train of true in combine_flag should be combined

  % start a block of flags to combine
  find(diff(combine_flag) > 0)
  
  find(diff(combine_flag) < 0)

  if (combine_flag(1)) 
  diff(combine_flag) 

  % remove events that are less than min_duration long

end

% time of center of energy
% rms
% mean frequency
 
