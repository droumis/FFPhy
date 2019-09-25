function bursts = measure_bursts(timestamp,timerange,threshold);
%MEASURE_BURSTS Identify bursts of a point process.
%
%   BURSTS = MEASURE_BURSTS(TIMESTAMP,TIMERANGE,THRESHOLD) computes various
%   measures related to the inter-event intervals of the point process events
%   with times TIMESTAMP observed during the time intervals given by TIMERANGE.
%
%   TIMESTAMP is a vector of real uint32 timestamps.
%
%   TIMERANGE is a set of uint32 time intervals (of the type validated by
%   IS_TIMERANGE). All event times in TIMESTAMP must lie within TIMERANGE.
%
%   THRESHOLD is a positive real scalar parameter; successive events which occur
%   less than THRESHOLD seconds apart are grouped together as bursts. A value of
%   6 ms is hisorically used to detect complex spikes of principal neurons in
%   hippocampus (Ranck, 1973). Larger values (e.g. 15 ms) detect "pseudo-bursts"
%   (Harris et al., 2001).
%
%   BURSTS is a vector of the same size as TIMESTAMP in which each element is
%   labeled as 0, 1, 2, 3, ... or NaN. 0 values indicate single events not
%   associated with a burst. 1, 2, 3, etc. values indicate the first, second,
%   third, etc. event within a burst. NaN indicates events that can not be
%   unambiguously determined (usually due to edge censoring at the start or end
%   of a time interval).
%
%   References:
%
%   Ranck J.B.J (1973) Studies on single units in dorsal hippocampal formation
%   and septum in unrestrained rats. I. Behavioral correlates and firing
%   repertoires. _Exp. Neurol._ 41:461-531.
%
%   Harris K., Hirase H., Leinekugel X., Henze D., Buzsáki G. (2001) Temporal
%   interaction between single spikes and complex spike bursts in hippocampal
%   pyramidal cells. _Neuron_ 32:141-149.
%
%   Cooper DC, Chung S, Spruston N. (2005) Output-mode transitions are
%   controlled by prolonged inactivation of sodium channels in pyramidal
%   neurons of subiculum. _PLoS Biology_ 3:e175.
%
%Depends on:
%   IS_TIMERANGE (written by SMK)
%   ISMEMBER_INTERVALS (written by SMK)
%
%Written by SMK, 2009 August 30.
%

TS_PER_SEC = 1e4;

if exist('is_timerange') ~= 2
  error('MEASURE_BURSTS depends on m-file IS_TIMERANGE (written by smk)');
end
if exist('ismember_intervals') ~= 2
  error('MEASURE_BURSTS depends on m-file ISMEMBER_INTERVALS (written by smk)');
end

if ~isnumeric(timestamp) || ~isvector(timestamp) || ...
    ~isa(timestamp,'uint32') || ~isreal(timestamp) || ...
    ~all(diff(timestamp) >= 0)
  error('TIMESTAMP must be a monotonically-increasing real uint32 vector');
end
if ~is_timerange(timerange)
  error(['TIMERANGE must be a valid set of non-overlapping uint32 time ' ...
      'intervals of the type validated by IS_TIMERANGE']);
end
if ~all(ismember_intervals(timestamp,timerange))
  error(['all elements of TIMESTAMP must lie within time intervals ' ...
      'specified by TIMERANGE']);
end
if ~isnumeric(threshold) || ~isscalar(threshold) || ~isreal(threshold) || ...
    (threshold <= 0) || isnan(threshold) || isinf(threshold)
  error('THRESHOLD must be a positive real scalar');
end
if (threshold > 1)
  warning(['THRESHOLD argument must be expressed in seconds; the value ' ...
      'that you specified (%f seconds) is unusually large. Are you sure ' ...
      'that this is correct?'],threshold);
end

% compute ISIs
diff_before = [NaN; diff(timestamp)];
diff_after = [diff(timestamp); NaN];
% correct at the edges of each time interval when we don't know about the
% preceding or following censored events and thus can not reliably assign an ISI
for i = 1:size(timerange,1)
  j_first = find((timestamp >= timerange(i,1)) & (timestamp < timerange(i,2)), ...
      1,'first');
  if (timestamp(j_first) - timerange(i,1) >= TS_PER_SEC*threshold)
    % if we observe a silence longer than threshold from the start of the time
    % interval to the first event within that interval, then we can be confident
    % that the preceding ISI is longer than the threshold. Inf is used here as a
    % sentinel value
    diff_before(j_first) = Inf;
  else
    % otherwise, we do not know whether the immediately preceding event was
    % within threshold
    diff_before(j_first) = NaN;
  end
  j_last = find((timestamp >= timerange(i,1)) & (timestamp < timerange(i,2)), ...
      1,'last');
  % the same logic applies for the last event within the time interval
  if (timerange(i,2) - timestamp(j_last) >= TS_PER_SEC*threshold)
    diff_after(j_last) = Inf;
  else
    diff_after(j_last) = NaN;
  end
end

bursts = nan(size(timestamp));
% single events are labeled with 0; note that NaN fails the
% greater-than-or-equal comparison
bursts((diff_before >= TS_PER_SEC*threshold) & ...
    (diff_after >= TS_PER_SEC*threshold)) = 0;
% the first event in each burst is labeled with 1; note that NaN fails the
% less-than comparison and also fails the greater-than-or-equal comparison
bursts((diff_before >= TS_PER_SEC*threshold) & ...
    (diff_after < TS_PER_SEC*threshold)) = 1;
% to fill in the remaining values, scan from the first event of each burst
burst_idx = find(bursts == 1);
for i = 1:numel(burst_idx)
  j = burst_idx(i);
  while true
    j = j + 1;
    % note that NaN always fails the less-than comparison
    % TODO: are there any special cases that this algorithm misses?!
    if (diff_before(j) < TS_PER_SEC*threshold)
      bursts(j) = bursts(j-1) + 1;
    else
      break;
    end
  end
end



