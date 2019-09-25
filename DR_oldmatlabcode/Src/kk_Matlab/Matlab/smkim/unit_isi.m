function [isi_pre, isi_post, burstflag] = unit_isi(unit,threshold)
%UNIT_ISI Compute inter-spike interval distribution and identify bursts.
%
%   [ISI_PRE, ISI_POST, BURSTFLAG] = UNIT_ISI(UNIT,THRESHOLD) takes single-unit
%   spike data UNIT and a scalar parameter THRESHOLD and returns three vectors.
%   UNIT must be a scalar struct of valid single-unit spike data of the type
%   validated by IS_UNIT, whose 'timestamp' field is a column vector of strictly
%   monotonically-increasing uint32 timestamps. THRESHOLD must be a real
%   floating-point scalar greater than zero. ISI_PRE, ISI_POST, BURSTFLAG are
%   real double column vectors of the same size as UNIT.timestamp. ISI_PRE
%   contains time lags to the preceding spike (in seconds); the first element
%   is NaN. ISI_POST contains positive time lags to the following spike (in
%   seconds); the last element is NaN. BURSTFLAG contains integers, where 0
%   means that the spike is a single spike, 1 means that the spike is the first
%   in a burst, 2 means that the spike is the second in a burst, etc. NaN
%   values are inserted where the assignment of bursts is ambiguous due to
%   cutoff at the start or end of UNIT.timerange.
%
%   THRESHOLD is a double scalar which corresponds to the inter-spike interval
%   that defines bursts, expressed in seconds. Successive events which occur
%   less than THRESHOLD apart are grouped together as bursts. A value of 6 ms
%   (=0.0060 seconds) is hisorically used to detect complex spikes of principal
%   neurons in hippocampus (Ranck, 1973). Larger values (e.g. 0.0150 s) capture
%   "pseudo-bursts" (Harris et al., 2001).
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
%   IS_UNIT (written by SMK)
%   ISMEMBER_INTERVALS (written by SMK)
%
%Written by SMK, 2009 August 30.
%

TS_PER_SEC = 1e4;

if exist('is_unit') ~= 2
  error('UNIT_ISI depends on m-file IS_UNIT (written by smk)');
end
if exist('ismember_intervals') ~= 2
  error('UNIT_ISI depends on m-file ISMEMBER_INTERVALS (written by smk)');
end

if ~is_unit(unit) || ~isscalar(unit)
  error('UNIT must be a scalar struct of single-unit spike data');
end
if ~isfloat(threshold) || ~isscalar(threshold) || ~isreal(threshold) || ...
    (threshold <= 0) || isnan(threshold) || isinf(threshold)
  error('THRESHOLD must be a positive real double scalar');
end
if (threshold > 1)
  warning(['THRESHOLD argument must be expressed in seconds; the value ' ...
      'that you specified (%f seconds) is unusually large. Are you sure ' ...
      'that this is correct?'],threshold);
end

% convert timestamps to seconds
t = double(unit.timestamp)/TS_PER_SEC;

if isempty(t)
  isi_pre = zeros([0,1]);
  isi_post = zeros([0,1]);
  burstflag = zeros([0,1]);
  return;
end

isi_pre = [NaN; diff(t)];
isi_post = [diff(t); NaN];
burstflag = nan(size(t));

% Single events are labeled with 0. Notice that NaN fails the
% greater-than-or-equal comparison
burstflag((isi_pre >= threshold) & (isi_post >= threshold)) = 0;

% The first event in each burst is labeled with 1.
burstflag((isi_pre >= threshold) & (isi_post < threshold)) = 1;

% Handle the start as best as possible
if (t(1) - double(unit.timerange(1))/TS_PER_SEC >= threshold)
  if (isi_post(1) >= threshold)
    burstflag(1) = 0;
  else
    burstflag(1) = 1;
  end
end

% Label the remaining spikes in each burst
burst_idx = find(burstflag == 1);
for i = 1:numel(burst_idx)
  j = burst_idx(i);
  assert(burstflag(j) == 1);
  while (j < numel(burstflag))
    if (isi_post(j) < threshold)
      burstflag(j+1) = burstflag(j) + 1;
      j = j + 1;
    else
      break;
    end
  end
end

% Now handle the end as best as possible
if isnan(burstflag(end)) && ...
    (double(unit.timerange(end))/TS_PER_SEC - t(end) >= threshold)
  burstflag(end) = 0;
end


