function [channel_config, datatype] = read_channel_config(filename)
%READ_CHANNEL_CONFIG Read recording channel info from NSpike config file.
%
%   [CHANNEL_CONFIG, DATATYPE] = READ_CHANNEL_CONFIG(CONFIG_FILENAME) reads the
%   NSpike config file specified by CONFIG_FILENAME and returns a
%   CHANNEL_CONFIG struct array which summarizes the recording configuration of
%   every channel listed in the config file and a DATATYPE string which
%   describes the type of data that is recorded by these channels. 
%
%   CONFIG_FILENAME must include the full filesystem path. Also, an error will
%   be raised if the config file does not contain a declaration of either
%   'datatype[ \t]+ SPIKE' or 'datatype[ \t]+ CONTINUOUS', or if the config
%   file seems to be invalid/inconsistent.
%
%   Each struct in CHANNEL_CONFIG corresponds to one recording channel listed in
%   the config file (presumably corresponding to one subject on a single day).
%   The fields are the following: 
%     electrode: index of the electrode
%     channel: index of the channel within the electrode group
%     depth: penetration depth of the tetrode in NSpike units of 1/12 turn
%     reference: 1x2 vector specifying reference channel [electrode channel]
%     passband: 1x2 vector specifying filter bandpass [low high] (in Hz)
%     threshold: spike-detection threshold (in microvolts)
%     dspnum: index of the NSpike DSP that processed the channel
%     dspchan: index of the NSpike DSP channel that was mapped to this channel
%     Fs: voltage sampling rate
%     source: CONFIG_FILENAME
%
%   DATATYPE is a string, either 'SPIKE' or 'CONTINUOUS', which describes the
%   type of the recording channels.
%
%   See also READ_SPIKE_CONFIG, READ_CONTINUOUS_CONFIG written by smk.
%
% Written by smk, 2009 February 05.
%

% Verify that filename specifies a real file and includes path
if ~ischar(filename)
  error('filename must be a string');
elseif (exist(filename) ~= 2)
  error('filename %s does not refer to a valid file on search path',filename);
elseif ~isdir(fileparts(filename))
  error('filename %s does not include path',filename);
elseif isempty(regexp(filename,'^.+(?=\.config$)','match','once'))
  error('file %s does not have *.config suffix',filename);
end
d = dir(filename);
if isempty(d)
  error('no files found that match file specifier %s',filename);
end
if (length(d) > 1)
  error('more than one file matches file specifier %s; be more specific', ...
      filename);
end
% guard against memory crash
FILE_TOO_BIG = 1e5; 
if (d.bytes > FILE_TOO_BIG)
  error(['file %s is too big to be a config file; please provide ' ...
      'correct config file names'],filename);
end
fid = fopen(filename,'r');
if (fid == -1)
  error('Could not open file %s',filename);
end
try
  % Note: the start-of-header sequence ought to be '%%BEGINCONFIG', but the gzip
  % compression sometimes corrupts the first character
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
    error('File %s has more than one "datatype SPIKE" declaration', ...
        filename);
  end
  if (numel(continuous_matches) > 1)
    error('File %s has more than one "datatype CONTINUOUS" declaration', ...
        filename);
  end
  if (~isempty(spike_matches) && ~isempty(continuous_matches))
    error('File %s specifies both SPIKE and CONTINUOUS datatypes',filename);
  elseif ~isempty(spike_matches)
    datatype = 'SPIKE';
  elseif ~isempty(continuous_matches)
    datatype = 'CONTINUOUS';
  else
    error(['File %s contains neither SPIKE nor CONTINUOUS datatype ' ...
        'declaration'],filename);
  end
catch
  error('Could not read file %s: %s',filename,ferror(fid));
end
if (fclose(fid) == -1)
  error('could not close file %s',filename);
end

% keep track of the sampling rate of each DSP, specified in lines of the form
% "dspinfo [dspnum] [samprate] [possibly some extra cruft]"
dspinfo = struct('dspnum',{},'samprate',{});
try
  % tmp_dspinfo is sorted by dspnum (first column)
  tmp_dspinfo = sortrows([str2double(regexp(txt, ...
      '(?<=\sdspinfo[ \t]+)\d+(?=[ \t]+\d+)','match')'), ...
      str2double(regexp(txt, ...
      '(?<=\sdspinfo[ \t]+\d+[ \t]+)\d+(?=[ \t]+)','match')')],1);
  assert(~isempty(tmp_dspinfo));
catch
  error('Config file %s does not contain expected dspinfo lines',filename);
end
num_dsps = size(tmp_dspinfo,1);
if (num_dsps ~= size(unique(tmp_dspinfo,'rows'),1))
  error('Found duplicate dspinfo lines for the same dsp in %s',filename);
end
for k = 1:num_dsps
  dspinfo(end+1,1) = struct('dspnum',tmp_dspinfo(k,1),'samprate', ...
      tmp_dspinfo(k,2));
end

% parse the individual channel declarations
try
  % dspnum is declared as "dspnum %d"
  dspnums = str2double(regexp(txt, ...
      '(?<=\sdspnum[ \t]+)\d+(?=[ \t]*\n)','match'))';
  % look up the corresponding entry in dspinfo to get sampling rate
  for k = 1:numel(dspnums)
    dspinfo_idx = find([dspinfo(:).dspnum] == dspnums(k));
    if (numel(dspinfo_idx) ~= 1)
      error('File %s does not have dspinfo for dspnum %d',filename,dspnums(k));
    end
    samprates(k,1) = dspinfo(dspinfo_idx).samprate;
  end
  % dspchan is declared as "dspchan %d"
  dspchans = str2double(regexp(txt, ...
      '(?<=\sdspchan[ \t]+)\d+(?=[ \t]*\n)','match'))';
  % electrode group number is declared as "number %d"
  electrodes = str2double(regexp(txt, ...
      '(?<=\snumber[ \t]+)\d+(?=[ \t]*\n)','match'))';
  % channel numbers within each group is declared as "electchan %d"
  channels = str2double(regexp(txt, ...
      '(?<=\selectchan[ \t]+)\d+(?=[ \t]*\n)','match'))';
  % depth is declared as "depth %d"
  depths = str2double(regexp(txt, ...
      '(?<=\sdepth[ \t]+)\d+(?=[ \t]*\n)','match'))';
  % reference electrode is declared as "refelect %d" and "refchan %d"
  references = [ ...
      str2double(regexp(txt, ...
      '(?<=\srefelect[ \t]+)\d+(?=[ \t]*\n)','match'))' , ...
      str2double(regexp(txt, ...
      '(?<=\srefchan[ \t]+)\d+(?=[ \t]*\n)','match'))' ];
  % threshold is declared as "thresh %d"
  thresholds = str2double(regexp(txt, ...
      '(?<=\sthresh[ \t]+)\d+(?=[ \t]*\n)','match'))';
  % filter passband is declared as "filter %d %d" (I could parse both %d in
  % the same regex, but it's onerous to repackage the resulting cell array
  % of strings into a matrix)
  passbands = [ ...
      str2double(regexp(txt, ...
      '(?<=\sfilter[ \t]+)\d+(?=[ \t]+\d+[ \t]*\n)','match'))' , ...
      str2double(regexp(txt, ...
      '(?<=\sfilter[ \t]+\d+[ \t]+)\d+(?=[ \t]*\n)','match'))' ];
catch
  error('Malformed config file %s',filename);
end
if size(dspchans,1) ~= size(samprates,1) || ...
    size(dspchans,1) ~= size(dspnums,1) || ...
    size(dspchans,1) ~= size(electrodes,1) || ...
    size(dspchans,1) ~= size(channels,1) || ...
    size(dspchans,1) ~= size(depths,1) || ...
    size(dspchans,1) ~= size(references,1) || ...
    size(dspchans,1) ~= size(thresholds,1) || ...
    size(dspchans,1) ~= size(passbands,1)
  error('Malformed config file %s',filename);
end

% populate channel_config. is there a less kludgy way to accomplish this?
channel_config = struct( ...
    'electrode'          , {}, ...
    'channel'            , {}, ...
    'depth'              , {}, ...
    'reference'          , {}, ...
    'passband'           , {}, ...
    'threshold'          , {}, ...
    'dspchan'            , {}, ...
    'dspnum'             , {}, ...
    'Fs'                 , {}, ...
    'source'             , {}); 
for k = 1:numel(dspchans)
  channel_config(end+1,1) = struct( ...
      'electrode'          , electrodes(k)  , ...
      'channel'            , channels(k)    , ...
      'depth'              , depths(k)      , ...
      'reference'          , references(k,:), ...
      'passband'           , passbands(k,:) , ...
      'threshold'          , thresholds(k)  , ...
      'dspchan'            , dspchans(k)    , ...
      'dspnum'             , dspnums(k)     , ...
      'Fs'                 , samprates(k)   , ...
      'source'             , filename       );
end
if isempty(channel_config)
  error('could not find any channel configuration details');
end

% check for duplicate entries
if (numel(unique([channel_config(:).dspchan])) ~= numel(channel_config))
  error('Found duplicate entries for the same dspchan in %s',filename);
end
if (size(unique([vertcat(channel_config(:).electrode), ...
    vertcat(channel_config(:).channel)],'rows'),1) ~= numel(channel_config))
  error('Found duplicate entries for the same electrode channel in %s', ...
      filename);
end
% For multiple channels on the same electrode, make sure that they agree on
% depth and reference
for i = 1:length(channel_config)
  same_electrode_idx = find([channel_config(:).electrode] == ...
      channel_config(i).electrode);
  depth = unique([channel_config(same_electrode_idx).depth]);
  if (numel(depth) ~= 1)
    error('Channels of the same electrode do not agree on depth in %s', ...
        filename);
  end
  reference = unique(vertcat(channel_config(same_electrode_idx).reference), ...
      'rows');
  if (size(reference,1) ~= 1)
    error('Channels of the same electrode do not share reference in %s', ...
        filename);
  end
end

