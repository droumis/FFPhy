function continuous = read_continuous(filename,continuous_config,session);
%READ_CONTINUOUS Read single-channel continuous data from an NSpike *.eeg file
%
%   CONTINUOUS = READ_CONTINUOUS(FILENAME,CONTINUOUS_CONFIG,SESSION) opens the
%   NSpike .eeg file specified by FILENAME and reads the stream of continuous
%   data that span the session(s) specified by SESSION, using the meta-data
%   given in CONTINUOUS_CONFIG.
%
%   SESSION must be a struct array with the following fields:
%         subject: a string
%             day: an integer
%           epoch: a string, descriptor of the session
%     environment: a string, descriptor of the environment
%          tstart: a time string of the form HH:MM:SS[.xxxx]
%            tend: a time string of the form HH:MM:SS[.xxxx]
%   The elements of SESSION must all have identical .subject and .day fields,
%   and the .day value must match the .day in CONTINUOUS_CONFIG.
%
%   CONTINUOUS_CONFIG must be a scalar struct with the following required fields:
%     subject: string, epoch of subject
%     day: integer, index of the recording day
%     electrode: index of the electrode
%     channel: index of the channel within the electrode group
%     region: string, epoch of the targeted brain region
%     depth: penetration depth of the electrode in NSpike units of 1/12 turn
%     reference: 1x2 vector specifying reference channel as [electrode channel]
%     passband: 1x2 array specifying filter bandpass in Hz
%
%   The return struct CONTINUOUS has the following fields:
%     subject: inherited from SESSION
%     day: inherited from SESSION/CONTINUOUS_CONFIG
%     epoch: inherited from SESSION
%     environment: inherited from SESSION
%     electrode: inherited from CONTINUOUS_CONFIG
%     channel: inherited from CONTINUOUS_CONFIG
%     depth: inherited from CONTINUOUS_CONFIG
%     hemisphere: inherited from CONTINUOUS_CONFIG
%     region: inherited from CONTINUOUS_CONFIG
%     reference: inherited from CONTINUOUS_CONFIG
%     passband: inherited from CONTINUOUS_CONFIG
%     Fs: the nominal sampling rate (double scalar)
%     timestamp: monotonically-increasing vector of uint32 timestamps
%     samples: a single-precision floating-point vector of voltages, expressed
%       in microvolts
%     sources: {FILENAME, CONTINUOUS_CONFIG.sources{:}}
%
%   See also IS_CONTINUOUS, written by smk.
%
%Depends on:
%   READ_CONTINUOUS_MEX (written by smk)
%   IS_CONTINUOUS_CONFIG (written by smk)
%   IS_SESSION (written by smk)
%   IS_CONTINUOUS (written by smk)
%
% Written by SMK, 2009 June 1.
%

if (exist('read_continuous_mex') ~= 3)
  error(['READ_CONTINUOUS depends on mex-file READ_CONTINUOUS_MEX ' ...
      '(written by smk)']);
end
if (exist('is_continuous_config') ~= 2)
  error(['READ_CONTINUOUS depends on m-file IS_CONTINUOUS_CONFIG ' ...
      '(written by smk)']);
end
if (exist('is_session') ~= 2)
  error('READ_CONTINUOUS depends on m-file IS_SESSION (written by smk)');
end
if (exist('is_continuous') ~= 2)
  error('READ_CONTINUOUS depends on m-file IS_CONTINUOUS (written by smk)');
end

% Compensate for unity-gain inversion by NSpike amplifiers
GAIN = -1;

% Verify that filename specifies a real file and includes path
if ~ischar(filename)
  error('filename must be a string');
elseif (exist(filename) ~= 2)
  error('filename %s does not refer to a valid file on search path',filename);
elseif ~isdir(fileparts(filename))
  error('filename %s does not include path',filename);
elseif isempty(regexp(filename,'^.+(?=\.eeg$)','match','once'))
  error('file %s does not have *.eeg suffix',filename);
end

% Check continuous_config struct *scalar*
if ~is_continuous_config(continuous_config) || ~isscalar(continuous_config)
  error(['CONTINUOUS_CONFIG does not appear to be a valid continuous config ' ...
      'struct scalar']);
end

% Check that SESSION is a valid session struct array
if ~is_session(session)
  error('SESSION does not appear to be a valid sessions data struct array');
end
% Check whether all elements of SESSION have matching subject and day
if (numel(unique({session(:).subject})) > 1)
  error('SESSION argument has elements whose subject fields differ');
end
if ~all([session(:).day] == continuous_config.day)
  error(['SESSION argument has elements whose day fields differ from the ' ...
      'day field of CONTINUOUS_CONFIG']);
end

% Initialize empty output, inheriting meta-data from SESSION and CONTINUOUS_CONFIG
continuous = rmfield(session,{'tstart','tend','timerange'});
[continuous(:).region]    = deal(continuous_config.region);
[continuous(:).hemisphere]= deal(continuous_config.hemisphere);
[continuous(:).depth]     = deal(continuous_config.depth);
[continuous(:).electrode] = deal(continuous_config.electrode);
[continuous(:).channel]   = deal(continuous_config.channel);
[continuous(:).reference] = deal(continuous_config.reference);
[continuous(:).passband]  = deal(continuous_config.passband);
[continuous(:).units]   = deal('microvolts');
[continuous(:).sources]   = deal({filename continuous_config.sources{:}});
  
for i = 1:numel(session)
  try
    tmp_continuous = read_continuous_mex(filename,session(i).timerange(1), ...
        session(i).timerange(2));
  catch
    error('read_continuous_mex failed for session %s',session(i).epoch);
  end
  if (numel(tmp_continuous) ~= 1)
    error('continuous file %s has one or more gaps between %s and %s', ...
        filename,session(i).tstart,session(i).tend);
  end
  % Check that the sampling rate returned by READ_CONTINUOUS_MEX matches the
  % sampling rate read from CONTINUOUS_CONFIG
  if (tmp_continuous.Fs ~= continuous_config.Fs)
    error(['continuous file %s suggests sampling rate %f, but ' ...
        'continuous_config gives sampling rate as %f'],filename, ...
        tmp_continuous.Fs,continuous_config.Fs);
  end
  % Divide samples by gain
  tmp_continuous.samples = tmp_continuous.samples / GAIN;
  % Add to output
  for j = 1:length(fieldnames(tmp_continuous))
    name = subsref(fieldnames(tmp_continuous),substruct('{}',{j}));
    continuous(i).(name) = tmp_continuous.(name);
  end
end

if ~is_continuous(continuous)
  error('There is a bug in either READ_CONTINUOUS or IS_CONTINUOUS');
end


