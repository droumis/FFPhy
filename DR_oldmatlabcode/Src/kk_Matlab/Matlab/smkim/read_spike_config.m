function spike_config = read_spike_config(config_filenames,recording_sites)
%READ_SPIKE_CONFIG Read spike-recording tetrode info from NSpike config file(s).
%
%   SPIKE_CONFIG = READ_SPIKE_CONFIG(FILENAMES,RECORDING_SITES) reads the
%   NSpike config files specified by FILENAMES and the RECORDING_SITES struct
%   array, and returns a SPIKE_CONFIG struct array which summarizes the
%   recording configuration of tetrodes. FILENAMES can be either a string or a
%   cell array of strings, which may contain directory paths and Unix-style
%   globs. Config files that do not contain 'datatype[ \t]+ SPIKE' will be
%   ignored with warning. Errors will be raised if any files that match
%   FILENAME are found that do not seem to be valid config files.  
%
%   RECORDING_SITES must be a struct array with the following fields:
%     subject: string naming the subject; must be identical for all elements
%     day: integer for the recording day; must be identical for all elements
%     electrode: integer index of the electrode
%     hemisphere: string description of the hemisphere (e.g. 'left','right','NA');
%     region: string description of the recording site  
%   An error will be raised if RECORDING_SITES does not contain an element
%   whose .electrode field matches each electrode listed in the config file.
%
%   Each struct in SPIKE_CONFIG corresponds to one tetrode in one subject on a
%   single day. The fields are the following:
%     subject: inherited from RECORDING_SITES
%     day: inherited from RECORDING_SITES
%     tetrode: index of the tetrode
%     depth: penetration depth of the tetrode in NSpike units of 1/12 turn
%     hemisphere: inherited from RECORDING_SITES
%     region: inherited from RECORDING_SITES
%     reference: 1x2 vector specifying reference channel as [electrode channel]
%     passbands: 4x2 array specifying filter bandpass in Hz for each channel;
%       each row is [low high]
%     thresholds: 4x1 vector of spike-detection thresholds (in microvolts)
%     dspnums: 4x1 vector of indices of the NSpike DSPs that processed the
%       channels
%     dspchans: 4x1 vector of indices of the NSpike DSP channels that were
%       mapped to the tetrode channels
%     Fs: sampling rate
%     sources: 2-element cell array of strings, containing names of the original
%       files that were sourced to generate this structure
%
%Depends on:
%   IS_RECORDING_SITES (written by smk)
%   READ_CHANNEL_CONFIG (written by smk)
%   IS_SPIKE_CONFIG (written by smk)
%
%Written by smk, 2009 February 05.
%

if (exist('is_recording_sites') ~= 2)
  error(['READ_SPIKE_CONFIG depends on m-file IS_RECORDING_SITES ' ...
      '(written by smk)']);
end
if (exist('read_channel_config') ~= 2)
  error(['READ_SPIKE_CONFIG depends on m-file READ_CHANNEL_CONFIG ' ...
      '(written by smk)']);
end
if (exist('is_spike_config') ~= 2)
  error(['READ_SPIKE_CONFIG depends on m-file IS_SPIKE_CONFIG ' ...
      '(written by smk)']);
end

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
    if (numel(spike_matches) > 1)
      error('config file %s has more than one "datatype SPIKE" declaration', ...
          fname);
    end
    if isempty(spike_matches)
      error('config file %s does not contain "datatype SPIKE" declaration', ...
          fname);
    end
    if (~isempty(spike_matches) && ~isempty(continuous_matches))
      error('config file %s specifies both SPIKE and CONTINUOUS datatypes', ...
          fname);
    end
    frewind(fid);
  catch
    error('Could not find SPIKE datatype in file %s: %s',fname,ferror(fid));
  end
  if (fclose(fid) == -1)
    error('could not close file %s',fname);
  end
  try
    [tmp_channel_config, datatype] = read_channel_config(fname);
  catch
    error('read_channel_config(%s) failed',fname);
  end
  if strcmp(datatype,'SPIKE')
    channel_config{i} = tmp_channel_config;
  else
    error('config file %s does not contain SPIKE datatype',fname);
  end
end
channel_config = vertcat(channel_config{:});
% check whether channel_config was populated
if isempty(channel_config)
  error('No channel configuration was found');
end

spike_config = struct( ...
    'subject'           , {}, ...
    'day'               , {}, ...
    'tetrode'           , {}, ...
    'depth'             , {}, ...
    'hemisphere'        , {}, ...
    'region'            , {}, ...
    'reference'         , {}, ...
    'passbands'         , {}, ...
    'thresholds'        , {}, ...
    'dspnums'           , {}, ...
    'dspchans'          , {}, ...
    'Fs'                , {}, ...
    'sources'           , {});
tetrode_nums = unique([channel_config(:).electrode]);
for i = 1:numel(tetrode_nums)
  channel_idx = find([channel_config(:).electrode] == tetrode_nums(i));
  % Check that four channels, indexed 0,1,2,3, are included in the tetrode
  [channels, channel_sort_order] = sort([channel_config(channel_idx).channel]);
  % Sort channel_idx so that the order corresponds to the channel order
  channel_idx = channel_idx(channel_sort_order);
  if ~(numel(channels) == 4) || ~all(channels == 0:3)
    error('The four channels of the tetrode are not correctly specified');
  end
  % Check that all four channels agree on depth, reference and samprate
  depth = unique([channel_config(channel_idx).depth]);
  if ~isscalar(depth)
    error('The four channels of the tetrode do not agree on depth');
  end
  reference = unique(vertcat(channel_config(channel_idx).reference),'rows');
  if (size(reference,1) ~= 1)
    error('The four channels of the tetrode do not agree on reference');
  end
  samprate = unique([channel_config(channel_idx).Fs]);
  if (size(samprate,1) ~= 1)
    error('The four channels of the tetrode do not agree on sampling rate');
  end
  % Read passband, threshold, dpsnum, dspchan for the four channels
  thresholds = vertcat(channel_config(channel_idx).threshold);
  passbands = vertcat(channel_config(channel_idx).passband);
  if any(passbands(:,end) <= passbands(:,1))
    error('Filter pass bands are specified in wrong order');
  end  
  dspnums = vertcat(channel_config(channel_idx).dspnum);
  dspchans = vertcat(channel_config(channel_idx).dspchan);
  % Look up corresponding entry in recording_sites to get .region string
  lookup_idx = find([recording_sites(:).electrode] == tetrode_nums(i));
  if isempty(lookup_idx)
    error('RECORDING_SITES has no matching entry for electrode %d', ...
        tetrode_nums(i));
  elseif ~isscalar(lookup_idx)
    error('more than one element of RECORDING SITES matches the same tetrode');
  end
  hemisphere = recording_sites(lookup_idx).hemisphere;
  region = recording_sites(lookup_idx).region;
  % Look up the filenames of the NSpike config file and the flat text record of
  % recording sites
  channel_config_source = unique({channel_config(channel_idx).source});
  % (remember that UNIQUE returns a cell array)
  if numel(channel_config_source) > 1
    error('The four channels of the tetrode do not agree on config source');
  end
  sources = {channel_config_source{1}, recording_sites(lookup_idx).source};
  spike_config(end+1,1) = struct( ...
      'subject'           , subject         , ...
      'day'               , day             , ...
      'tetrode'           , tetrode_nums(i) , ...
      'depth'             , depth           , ...
      'hemisphere'        , hemisphere      , ...
      'region'            , region          , ...
      'reference'         , reference       , ...
      'passbands'         , passbands       , ...
      'thresholds'        , thresholds      , ...
      'dspnums'           , dspnums         , ...
      'dspchans'          , dspchans        , ...
      'Fs'                , samprate        , ...
      'sources'           , {sources}       );
      % we have to wrap sources in curly braces because of a stupid "feature" of
      % the STRUCT function
end  
 
if ~is_spike_config(spike_config)
  error('There is a bug in either READ_SPIKE_CONFIG or IS_SPIKE_CONFIG');
end

  
