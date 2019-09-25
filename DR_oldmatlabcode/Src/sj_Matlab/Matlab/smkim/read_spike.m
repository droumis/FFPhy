function [spike, num_rejected] = read_spike(filename,spike_config,amplitude_cutoff)
%READ_SPIKE Readd tetrode spike data from NSpike *.tt packed binary format.
%
%   SPIKE = READ_SPIKE(FILENAME,SPIKE_CONFIG,AMPLITUDE_CUTOFF) reads the
%   treshold-triggered snippets of waveform samples from the packed binary
%   data format in FILENAME and returns a SPIKE struct. 
%
%   FILENAME must specify an NSpike *.tt file of the type outputted by
%   nspike_extract -spike and must include the path to the file.
%
%   SPIKE_CONFIG must be a scalar struct with the following required fields:
%     subject: string, name of subject
%     day: integer, index of the recording day
%     tetrode: index of the tetrode
%     region: string, name of the targeted brain region
%     depth: penetration depth of the tetrode in NSpike units of 1/12 turn
%     reference: 1x2 vector specifying reference channel as [electrode channel]
%     passbands: 4x2 array specifying filter bandpass in Hz for each channel;
%       each row is [low high]
%     thresholds: 4x1 vector of spike-detection thresholds (in microvolts)
%     Fs: sampling rate
%   An error will be raised if any of the spikes does not exceed its
%   respective threshold at any sample in the window.
%
%   AMPLITUDE_CUTOFF must be a 4-element vector of positive real values. A
%   threshold-trigger event will be censored if the waveform recorded on *any*
%   of the channels 1-4 exceed its respective cutoff in AMPLITUDE_CUTOFF(1:4).
%   Note that this behavior can result in the accidental rejection of valid
%   spikes if one of the channels is noisy and its AMPLITUDE_CUTOFF is not set
%   to Inf. An appropriate AMPLITUDE_CUTOFF can be set interactively using the
%   function CHECK_TT.
%
%   SPIKE is struct with the following required fields:
%     subject: inherited from SPIKE_CONFIG
%     day: inherited from SPIKE_CONFIG
%     tetrode: inherited from SPIKE_CONFIG
%     depth: inherited from SPIKE_CONFIG
%     hemisphere: inherited from SPIKE_CONFIG
%     region: inherited from SPIKE_CONFIG
%     reference: inherited from SPIKE_CONFIG
%     passbands: inherited from SPIKE_CONFIG
%     thresholds: inherited from SPIKE_CONFIG
%     timestamp: a Nx1 vector of uint32 (in units of 1e-4 seconds)
%     samples: 40x4xN 3-dimensional array of int16; samples(:,I,J) is the
%       vector of 40 samples recorded on the Ith channel during the Jth
%       threshold-trigger event
%     Fs: inherited from SPIKE_CONFIG
%     samples_before_trigger: number of samples that were saved before the
%       threshold-crossing; equals 8 for an NSpike *.tt file
%     amplitude_cutoff: AMPLITUDE_CUTOFF
%     sources: {FILENAME, SPIKE_CONFIG.sources{:}}
%
%   [SPIKE, NUM_REJECTED] = READ_SPIKE(FILENAME,THRESHOLDS,AMPLITUDE_CUTOFF)
%   reports the number of events that were excluded for violating
%   AMPLITUDE_CUTOFF. 
%
%   See also IS_SPIKE and CHECK_TT, written by smk.
%
%Depends on:
%   READ_BINARY_RECORDS (written by smk)
%   IS_SPIKE_CONFIG (written by smk)
%   IS_SPIKE (written by smk)
%
%Written by smk, 2009 February 13
%

if (exist('read_binary_records') ~= 2)
  error('READ_SPIKE depends on m-file READ_BINARY_RECORDS (written by smk)');
end
if (exist('is_spike_config') ~= 2)
  error('READ_SPIKE depends on m-file IS_SPIKE_CONFIG (written by smk)');
end
if (exist('is_spike') ~= 2)
  error('READ_SPIKE depends on m-file IS_SPIKE (written by smk)');
end

% Number of channels per tetrode
NCHANNELS = 4;
% Number of samples per window in the original tt file
NPRE_THRESH_POINTS = 8;
NPOST_THRESH_POINTS = 32;
NPOINTS_PER_SPIKE = NPRE_THRESH_POINTS + NPOST_THRESH_POINTS;
% NSpike acquires voltage samples with (-1) gain due to inverting op-amps
GAIN = int16(-1);

if ~isnumeric(amplitude_cutoff) || ~isvector(amplitude_cutoff) || ...
    ~isreal(amplitude_cutoff) || ~all(amplitude_cutoff >= 0) || ...
    (numel(amplitude_cutoff) ~= NCHANNELS)
  error('amplitude_cutoff must be a %d-element vector of positive reals', ...
      NCHANNELS);
elseif (size(amplitude_cutoff,1) > 1)
  % Reshape as a row vector
  amplitude_cutoff = amplitude_cutoff';
end

% Verify that filename specifies a real file and includes path
if ~ischar(filename)
  error('filename must be a string');
elseif (exist(filename) ~= 2)
  error('filename %s does not refer to a valid file on search path',filename);
elseif ~isdir(fileparts(filename))
  error('filename %s does not include path',filename);
elseif isempty(regexp(filename,'^.+(?=\.tt$)','match','once'))
  error('file %s does not have *.tt suffix',filename);
end

% Check that spike_config is a valid struct *scalar*
if ~is_spike_config(spike_config) || ~isscalar(spike_config)
  error(['SPIKE_CONFIG does not appear to be a valid ' ...
      'spike config data struct scalar']);
end

% Initialize and inherit meta-data
spike = rmfield(spike_config,{'dspnums','dspchans'});
% Convert thresholds to row vector
if (size(spike.thresholds,2) < size(spike.thresholds,1))
  spike.thresholds = spike.thresholds';
end
spike.samples_before_trigger = NPRE_THRESH_POINTS;
spike.amplitude_cutoff = amplitude_cutoff;
spike.timestamp = uint32(zeros([0 1]));
spike.samples = int16(zeros([NCHANNELS NPOINTS_PER_SPIKE 0]));
spike.sources = {filename spike_config.sources{:}};

% read data from tt file using READ_BINARY_RECORDS
recformat = cell2struct({ ...
    'timestamp', 'uint32', 1; ...
    'samples'  , 'int16' , NCHANNELS*NPOINTS_PER_SPIKE }, ...
    {'name','type','count'},2);
eoh = '%%ENDHEADER\n';
try
  tmp = read_binary_records(filename,eoh,recformat,Inf);
  spike.timestamp = tmp.timestamp;
  num_events = numel(spike.timestamp);
  % samples field has num_events rows and (NCHANNELS*NPOINTS_PER_SPIKE) columns;
  % each row contains samples in interlaced order, e.g. 
  %
  % (channel1,sample1), (channel2,sample1), ... (channelNCHANNELS,sample1), 
  % (channel1,sample2), (channel2,sample2), ... (channelNCHANNELS,sample2), 
  % (channel1,sample3), ...
  % (channel1,sampleNPOINTS_PER_SPIKE), (channel2,sampleNPOINTS_PER_SPIKE), ...
  %
  % We want to convert this to a 3-dimensional array of size 
  % NPOINTS_PER_SPIKE * NCHANNELS * num_events. Also, invert sign to undo the -1
  % gain in NSpike (so that, as expected, positive values mean positive recorded
  % extracellular voltage)
  spike.samples = permute(reshape(tmp.samples', ...
      [NCHANNELS NPOINTS_PER_SPIKE num_events]),[2 1 3]) / GAIN;
catch
  error('No valid data found in file %s',filename);
end

if (num_events == 0)
  return;
end

% Verify that thresholds are appropriate: in every event, at least one channel
% must exceed its threshold
below_threshold = find(~any(max(-spike.samples,[],1) >= ...
    repmat(spike.thresholds,[1 1 num_events])));
if ~isempty(below_threshold)
  warning(['%d (out of %d) events in %s fail to exceed thresholds %s. ' ...
      'THESE WILL BE CENSORED!'],numel(below_threshold),num_events, ...
      filename,mat2str(spike.thresholds));
  spike.samples(:,:,below_threshold) = [];
  spike.timestamp(below_threshold) = [];
  num_events = size(spike.samples,3);
end

% Remove events that exceed amplitude_cutoff
if all(amplitude_cutoff(:) == 0)
  spike.timestamp = uint32(zeros([0 1]));
  spike.samples = int16(zeros([NCHANNELS NPOINTS_PER_SPIKE 0]));
  return;
end
cull_idx = find(any(max(abs(spike.samples),[],1) > ...  
    repmat(amplitude_cutoff,[1 1 num_events]),2));
num_rejected = numel(cull_idx);
if ~isempty(cull_idx)
  disp(sprintf('%d events exceed amplitude cutoff %s, will be excluded', ...
      num_rejected,mat2str(amplitude_cutoff)));
  keep_idx = setdiff(1:num_events,cull_idx);
  spike.timestamp = spike.timestamp(keep_idx,1);
  spike.samples = spike.samples(:,:,keep_idx);
  if isempty(spike.timestamp)
    warning('spike struct contains no trigger events');
  end
end

if (numel(spike.timestamp) > 1)
  % Report diagnostic information to the user
  timediffs = diff(spike.timestamp);
  idx = find(timediffs <= 0);
  if ~isempty(idx)
    for i = 1:numel(idx)
      warning('Non-increasing timestamps between consecutive events: %d', ...
          timediffs(idx(i)));
    end
  end
end

if ~is_spike(spike)
  error('there is bug in either READ_SPIKE or IS_SPIKE');
end


