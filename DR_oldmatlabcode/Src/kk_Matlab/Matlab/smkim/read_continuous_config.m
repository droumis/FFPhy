function continuous_config = read_continuous_config(config_filenames,recording_sites)
%READ_CONTINUOUS_CONFIG Read continuous recording channel configuration from NSpike config file(s).
%
%   CONTINUOUS_CONFIG = READ_CONTINUOUS_CONFIG(FILENAMES,RECORDING_SITES) reads
%   the NSpike config files specified by FILENAMES and returns a
%   CONTINUOUS_CONFIG struct array which summarizes the configuration of
%   channels for recording continuous data (e.g. LFP, EMG). FILENAMES can be
%   either a string or a cell array of strings, which may contain directory
%   paths and Unix-style globs. Config files that do not contain 'datatype[ \t]+
%   CONTINUOUS' will be ignored with warning. Errors will be raised if any files
%   that match FILENAME are found that do not seem to be valid config files.
%
%   RECORDING_SITES must be a struct array with the following fields:
%     subject: string naming the subject; must be identical for all elements
%     day: integer for the recording day; must be identical for all elements
%     electrode: integer index of the electrode
%     hemisphere: string description of the hemisphere (e.g. 'left','right','NA')
%     region: string description of the recording site  
%   An error will be raised if RECORDING_SITES does not contain an element
%   whose .electrode field matches each electrode listed in the config file.
%
%   Each struct in CONTINUOUS_CONFIG corresponds to one continuous recording
%   channel in one subject on a single day. (It is assumed that the channel
%   settings are not altered throughout the day.) The fields are the following: 
%     subject: inherited from RECORDING_SITES
%     day: inherited from RECORDING_SITES
%     hemisphere: inherited from RECORDING_SITES
%     region: inherited from RECORDING_SITES
%     electrode: index of the electrode
%     channel: index of the channel within the electrode group
%     depth: penetration depth of the tetrode in NSpike units of 1/12 turn
%     reference: 1x2 vector specifying reference channel [electrode channel]
%     passband: 1x2 vector specifying filter bandpass [low high] (in Hz)
%     dspnum: index of the NSpike DSP that processed the channel
%     dspchan: index of the NSpike DSP channel that was mapped to this channel
%     Fs: voltage sampling rate
%     sources: 2-element cell array of strings, containing names of the original
%       files that were sourced to generate this structure
%
%   This function ignores the DSP channel that carries the video
%   synchronization square-wave (channel #63 in the standard NSpike setup).
%
%   See also IS_CONTINUOUS_CONFIG, written by smk.
%
%Depends on:
%   IS_RECORDING_SITES (written by smk)
%   READ_CHANNEL_CONFIG (written by smk)
%   IS_CONTINUOUS_CONFIG (written by smk)
%
%Written by smk, 2009 February 05.
%

if (exist('is_recording_sites') ~= 2)
  error(['READ_CONTINUOUS_CONFIG depends on m-file IS_RECORDING_SITES ' ...
      '(written by smk)']);
end
if (exist('read_channel_config') ~= 2)
  error(['READ_CONTINUOUS_CONFIG depends on m-file READ_CHANNEL_CONFIG ' ...
      '(written by smk)']);
end
if (exist('is_continuous_config') ~= 2)
  error(['READ_CONTINUOUS_CONFIG depends on m-file IS_CONTINUOUS_CONFIG ' ...
      '(written by smk)']);
end

% This is a NSpike hardware-specific constant
VIDEO_SYNC_CHANNEL = 63;

if ischar(config_filenames)
  config_filenames = {config_filenames};
end
if ~iscellstr(config_filenames)
  error('config_filenames must be a string or a cell array of strings');
end
if any(cellfun(@isempty,config_filenames))
  error('config_filenames can not be empty');
end

% Validate recording_sites argument
if ~is_recording_sites(recording_sites)
  error(['RECORDING_SITES does not appear to be a valid recording sites ' ...
      'data struct array']);
end
% Check whether recording_sites struct refers to exactly one subject/day; this
% condition is imposed because NSpike config files don't contain any information
% that distinguishes one day or subject from another
subject = unique({recording_sites(:).subject});
if (numel(subject) ~= 1)
  error('all elements of RECORDING_SITES must share same subject field');
else
  subject = subject{1};
end
day = unique([recording_sites(:).day]);
if ~isscalar(day)
  error('all elements of RECORDING_SITES must share same day field');
end

% Initialize empty cell array for storing channel config info read from
% multiple config files
channel_config = {};
for i = 1:length(config_filenames)
  fname = config_filenames{i};
  % Verify that filename specifies a real file and includes path
  if ~ischar(fname)
    error('fname must be a string');
  elseif (exist(fname) ~= 2)
    error('fname %s does not refer to a valid file on search path',fname);
  elseif ~isdir(fileparts(fname))
    error('fname %s does not include path',fname);
  elseif isempty(regexp(fname,'^.+(?=\.config$)','match','once'))
    error('file %s does not have *.config suffix',fname);
  end
  d = dir(fname);
  if isempty(d)
    error('no files found that match file specifier %s',fname);
  end
  if (length(d) > 1)
    error('more than one file matches file specifier %s; be more specific', ...
        fname);
  end
  % guard against memory crash
  FILE_TOO_BIG = 1e5; 
  if (d.bytes > FILE_TOO_BIG)
    error(['file %s is too big to be a config file; please provide ' ...
        'correct config file names'],fname);
  end
  fid = fopen(fname,'r');
  if (fid == -1)
    error('Could not open file %s',fname);
  end
  try
    % Note: the start-of-header sequence ought to be '%%BEGINCONFIG', but the
    % gzip compression sometimes corrupts the first character
    while isempty(regexp(fgets(fid),'%BEGINCONFIG\n'))
      if feof(fid)
        error('File %s does not have %BEGINCONFIG character sequence', ...
            filename);
      end
    end
    txt = fread(fid,Inf,'*char')';
    % look for datatype declaration
    continuous_matches = regexp(txt,'datatype[[ \t]+\w+]*[ \t]+CONTINUOUS');
    spike_matches = regexp(txt,'datatype[[ \t]+\w+]*[ \t]+SPIKE');
    if (numel(continuous_matches) > 1)
      error(['config file %s contains more than one "datatype CONTINUOUS" ' ...
          'declaration'],fname);
    end
    if isempty(continuous_matches)
      error(['config file %s does not contain a "datatype CONTINUOUS" ' ...
          'declaration'],fname);
    end
    if (~isempty(spike_matches) && ~isempty(continuous_matches))
      error('config file %s specifies both SPIKE and CONTINUOUS datatypes', ...
          fname);
    end
    frewind(fid);
  catch
    error('Could not read file %s: %s',fname,ferror(fid));
  end
  if (fclose(fid) == -1)
    error('could not close file %s',fname);
  end
  try
    [tmp_channel_config, datatype] = read_channel_config(fname);
  catch
    error('read_channel_config(%s) failed',fname);
  end
  if strcmp(datatype,'CONTINUOUS')
    channel_config{i} = tmp_channel_config;
  else
    error('config file %s does not contain CONTINUOUS datatype',fname);
  end
end
% check whether channel_config was populated
channel_config = vertcat(channel_config{:});
if isempty(channel_config)
  error('No channel configuration was found');
end
% remove the channels with the video sync signal (fake "LFP" channel)
videosync_idx = find([channel_config(:).dspchan] == VIDEO_SYNC_CHANNEL);
if ~isempty(videosync_idx)
  channel_config(videosync_idx) = [];
end

% copy channel_config to continuous_config, removing the threshold field
continuous_config = rmfield(channel_config,'threshold');
% Copy meta-data from RECORDING_SITES
[continuous_config.subject] = deal(subject);
[continuous_config.day] = deal(day);
for i = 1:numel(continuous_config)
  % Look up corresponding entry in recording_sites
  match_idx = find( ...
      [recording_sites(:).electrode] == continuous_config(i).electrode);
  if (numel(match_idx) ~= 1)
    error(['RECORDING_SITES has no matching entry for electrode %d, ' ...
        'channel %d'],continuous_config(i).electrode, ...
        continuous_config(i).channel);
  end
  continuous_config(i).hemisphere = recording_sites(match_idx).hemisphere;
  continuous_config(i).region = recording_sites(match_idx).region;
  % Combine the NSpike config filename with the filename of the flat text record
  % of recording sites
  continuous_config(i).sources = {continuous_config(i).source, ...
      recording_sites(match_idx).source};
end

% Remove the source (singular! not sourceS) field, because that information is
% now kept in the .sources field cell array
continuous_config = rmfield(continuous_config,'source');

if ~is_continuous_config(continuous_config)
  error(['There is a bug in either READ_CONTINUOUS_CONFIG or ' ...
      'IS_CONTINUOUS_CONFIG']);
end

  
