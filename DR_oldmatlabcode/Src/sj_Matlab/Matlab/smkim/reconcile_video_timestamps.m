function timestamps = reconcile_video_timestamps(cpudsptimecheck_filename,cpupostimestamp_filename,videosync_filename)
%RECONCILE_VIDEO_TIMESTAMPS Generate corrected video frame timestamps from NSpike data.
%
%   TIMESTAMPS = RECONCILE_VIDEO_TIMESTAMPS(CPUDSPTIMECHECK_FILENAME,
%   CPUPOSTIMESTAMP_FILENAME, VIDEOSYNC_FILENAME) reads NSpike timestamp outputs
%   to produce a vector of corrected timestamps for the video frames.
%
%   CPUDSPTIMECHECK_FILENAME must specify an NSpike CPU-DSP clock
%   synchronization reference file. This file contains a series of struct
%   records written as packed binary data. Each record contains a triplet of
%   uint32_t values: the local time of the computer (i.e. polling the
%   local cpu clock), then the (authoritative) time received when polling the
%   DSP, and lastly the local time of the cpu clock again.
%
%   CPUPOSTIMESTAMP_FILENAME must specify an NSpike cpu videoap timestamp file. 
%   This file contains a series of uint32_t values as packed binary data,
%   which are the timestamps of the video capture events as determined by
%   polling the local cpu clock.
%
%   VIDEOSYNC_FILENAME must specify an NSpike continuous-data file which
%   contains the video sync voltage sweep from the camera. The video sync signal
%   is a square wave. We want to detect the times when the voltage jumps down
%   from the (+) plateau to the (-) plateau. Because the signal is passed
%   through a low-pass causal filter, this transition is actually smoothed in
%   the recorded data. So really, we want to find the inflection points when the
%   signal goes from approximately flat at the (+) plateau to precipitously
%   decreasing. Note that this function will fail if the polarity of the video
%   sync signal is inverted or if the signal was not subject to a causal
%   low-pass filter.
%
%   The output TIMESTAMPS is a uint32 vector whose length is the number of video
%   capture events listed in the file specified by CPUPOSTIMESTAMP_FILENAME.
%
%Depends on:
%   READ_BINARY_RECORDS (written by smk)
%   READ_CONTINUOUS_MEX (written by smk)
%   ROBUSTFIT (MATLAB Statistics Toolbox)
%
%Written by smk, 2009 July 30.
%

% DO NOT TAMPER WITH THIS CONSTANT CONVERSION FACTOR!
TS_PER_SEC = 1e4;

% Maximum allowed polling delay between DSP time check before and DSP time check
% after; if this delay was exceeeded, then the computer timestamp for that frame
% is considered to be invalid and we rely on integrated video sync pulses
% instead.
MAX_CPU_DELAY = 10; % 1 millisecond

% Maximum allowed discrepancy between the time of the video sync signal
% transition and the expected time of the frame according to interpolation
MAX_ALLOWED_DISCREPANCY = 10; % 1 millisecond

if (exist('read_binary_records') ~= 2)
  error(['RECONCILE_VIDEO_TIMESTAMPS depends on m-file ' ...
      'READ_BINARY_RECORDS (written by smk)']);
end
if (exist('read_continuous_mex') ~= 3)
  error(['RECONCILE_VIDEO_TIMESTAMPS depends on mex-file ' ...
      'READ_CONTINUOUS_MEX (written by smk)']);
end
if (exist('robustfit') ~= 2)
  error(['RECONCILE_VIDEO_TIMESTAMPS depends on m-file ' ...
      'ROBUSTFIT (in MATLAB Statistics Toolbox)']);
end

% Verify that filenames are correct
if ~ischar(cpudsptimecheck_filename)
  error('cpudsptimecheck_filename must be a string');
elseif (exist(cpudsptimecheck_filename) ~= 2)
  error(['cpudsptimecheck_filename %s does not refer to a valid file on ' ...
      'search path'],cpudsptimecheck_filename);
elseif isempty(regexp(cpudsptimecheck_filename, ...
    '^.+(?=\.cpudsptimecheck$)','match','once'))
  error('file %s does not have *.cpudsptimecheck suffix', ...
      cpudsptimecheck_filename);
end
if ~ischar(cpupostimestamp_filename)
  error('cpupostimestamp_filename must be a string');
elseif (exist(cpupostimestamp_filename) ~= 2)
  error(['cpupostimestamp_filename %s does not refer to a valid file on ' ...
      'search path'],cpupostimestamp_filename);
elseif isempty(regexp(cpupostimestamp_filename, ...
    '^.+(?=\.cpupostimestamp$)','match','once'))
  error('file %s does not have *.cpupostimestamp suffix', ...
      cpupostimestamp_filename);
end
if ~ischar(videosync_filename)
  error('videosync_filename must be a string');
elseif (exist(videosync_filename) ~= 2)
  error(['videosync_filename %s does not refer to a valid file on ' ...
      'search path'],videosync_filename);
elseif isempty(regexp(videosync_filename, ...
    '^.+(?=\.videosync$)','match','once'))
  error('file %s does not have *.videosync suffix', ...
      videosync_filename);
end

try
  videosync = read_continuous_mex(videosync_filename); 
catch
  error('could not read video sync file %s');
end
% Identify transition points
for i = 1:length(videosync)
  t = double(videosync(i).timestamp);
  x = diff(videosync(i).samples) < 0;
  transitions{i} = t(intersect(find(x), 1 + find(~x)));
end
transitions = vertcat(transitions{:});

% cpudsptimecheck contains triplets of timestamps which were polled about once
% every 10 seconds during data acquisition. These are used for interpolating
% unreliable cpu timestamps with respect to the authoritative DSP timestamps
try
  cpudsptimecheck = read_binary_records(cpudsptimecheck_filename, ...
      '%%ENDHEADER\n',struct('name',{'cpu_before';'dsp_between'; ...
      'cpu_after'},'type',{'uint32';'uint32';'uint32'},'count',{1;1;1}),Inf);
catch
  error('Could not read cpu-dsp time checks in %s',cpudsptimecheck_filename);
end
if ( any(cpudsptimecheck.cpu_after < cpudsptimecheck.cpu_before) || ...
    any(diff(cpudsptimecheck.cpu_after) <= 0) || ...
    any(diff(cpudsptimecheck.cpu_before) <= 0) )
  error('consecutive DSP timestamps in %s are out of order.', ...
      cpudsptimecheck_filename);
end

% For each segment of continuous record in videosync, identify time check events
% for which the CPU and DSP times are implausibly discrepant (we expect them to
% exhibit a quasi-linear relationship corrupted by noise (not necessarily i.i.d.
% gaussian)). We can use ROBUSTFIT to detect outliers.
[junk, stats1] = robustfit(double(cpudsptimecheck.dsp_between), ...
    double(cpudsptimecheck.cpu_before),'bisquare');
[junk, stats2] = robustfit(double(cpudsptimecheck.dsp_between), ...
    double(cpudsptimecheck.cpu_after),'bisquare');
% Censor these anomalous time check events, as well as those which exceeded
% polling delay (indicative of a hiccup in the frame capture)
bad_timecheck_idx = find((stats1.w < 0.75) | (stats2.w < 0.75) | ...
    (cpudsptimecheck.cpu_after - cpudsptimecheck.cpu_before > MAX_CPU_DELAY));
if ~isempty(bad_timecheck_idx)
  warning('%d out of %d time checks in %s will be excluded', ...
      numel(bad_timecheck_idx),numel(cpudsptimecheck.dsp_between), ...
      cpudsptimecheck_filename);
  cpudsptimecheck.cpu_before(bad_timecheck_idx) = [];
  cpudsptimecheck.dsp_between(bad_timecheck_idx) = [];
  cpudsptimecheck.cpu_after(bad_timecheck_idx) = [];
end

% cpu_timestamps contains timestamps with reference to local computer clock,
% one per each video frame. The task is to translate these inaccurate computer
% timestamps into "correct" DSP-based timestamps
try
  cpu_timestamps = read_binary_records(cpupostimestamp_filename, ...
      '%%ENDHEADER\n',struct('name',{'cpu'},'type',{'uint32'}, ...
      'count',{1}),Inf);
  cpu_timestamps = cpu_timestamps.cpu; 
catch
  error('Could not read cpu timestamps in %s',cpupostimestamp_filename);
end

% The naive approach is linear interpolation using cpudsptimecheck, which
% roughly compensates for clock skew of compute relative to the DSP.
expected_timestamps = interp1(0.5*(double(cpudsptimecheck.cpu_before) + ...
    double(cpudsptimecheck.cpu_after)),double(cpudsptimecheck.dsp_between), ...
    double(cpu_timestamps),'linear','extrap');
% Look up the nearest video sync signal transition
timestamps = uint32(round(interp1(transitions,transitions, ...
    expected_timestamps,'nearest')));

% Check for gross discrepancies (usually, these result because the video sync
% waveform was clipped at the start or end of a data acquisition interval)
dt = double(timestamps) - expected_timestamps;
uncertain_frames = abs(dt) > MAX_ALLOWED_DISCREPANCY;
timestamps(uncertain_frames) = uint32(round( ...
    expected_timestamps(uncertain_frames)));
warning(['Interpolated timestamp was substituted for %d video frames out ' ...
    'of %d total. Median of discrepancies is %f ms. Estimated standard ' ...
    'deviation of discrepancies is %f ms'], ...
    nnz(uncertain_frames),numel(timestamps), ...
    median(double(timestamps) - expected_timestamps)/TS_PER_SEC*1e3, ...
    std(double(timestamps) - expected_timestamps)/TS_PER_SEC*1e3);
    
  
