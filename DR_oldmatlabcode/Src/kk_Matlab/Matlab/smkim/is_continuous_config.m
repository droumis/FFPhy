function bool = is_continuous_config(continuous_config)
%IS_CONTINUOUS_CONFIG Validate whether the input is an continguous_config struct array
%
%
%   IS_CONTINUOUS_CONFIG(CONTINUOUS_CONFIG) returns true if CONTINUOUS_CONFIG is
%   a struct array that contains valid continuous recording channel
%   configuration meta-data.
%   
%   See also READ_CONTINUOUS_CONFIG, written by smk.
%
%Written by SMK 2009 June 22.
%

% This is determined by the NSpike hardware.
DSPCHAN_ALLOWED = [
    0;  1;  2; ...
    4;  5;  6; ...
    8;  9; 10; ...
   12; 13; 14; ...
   16; 17; 18; ...
   20; 21; 22; ...
   24; 25; 26; ...
   28; 29; 30; ...
   32; 33; 34; ...
   36; 37; 38; ...
   40; 41; 42; ...
   44; 45; 46; ...
   48; 49; 50; ...
   52; 53; 54; ...
   56; 57; 58; ...
   60; 61; 62; ...
   64; 65; 66; ...
   68; 69; 70; ...
   72; 73; 74; ...
   (76:126)' ];

% This is determined by the NSpike hardware.
DSPNUM_ALLOWED = 1:8;

REQUIRED_FIELDS = { ...
    'subject'               , ...
    'day'                   , ...
    'electrode'             , ...
    'channel'               , ...
    'depth'                 , ...
    'hemisphere'            , ...
    'region'                , ...
    'reference'             , ...
    'passband'              , ...
    'dspnum'                , ...
    'dspchan'               , ...
    'Fs'                    , ...
    'sources'               };

if isempty(continuous_config) && isstruct(continuous_config) && ...
    all(isfield(continuous_config,REQUIRED_FIELDS))
  warning('CONTINUOUS_CONFIG struct is empty');
  bool = true;

elseif ~isstruct(continuous_config) || ...
    ~all(isfield(continuous_config,REQUIRED_FIELDS))
  warning('missing required field(s)');
  bool = false;

elseif ~iscellstr({continuous_config(:).subject}) || ...
    ~all(strcmp({continuous_config(:).subject}, ...
    regexp({continuous_config(:).subject},'[a-zA-Z0-9]*','match','once')))
  warning(['subject field must be a string containing characters ' ...
      '[a-zA-Z0-9]']);
  bool = false;

elseif ~all(cellfun(@isnumeric,{continuous_config(:).day})) || ...
    ~all(cellfun(@isscalar,{continuous_config(:).day})) || ...
    ~all(cellfun(@isreal,{continuous_config(:).day})) || ...
    ~all(cellfun(@isfinite,{continuous_config(:).day})) || ...
    ~all(cellfun(@(c) (c >= 0),{continuous_config(:).day})) || ...
    ~all(cellfun(@(c) (round(c) == c),{continuous_config(:).day}))
  warning('day field must be a real non-negative integer scalar');
  bool = false;

elseif ~all(cellfun(@isnumeric,{continuous_config(:).electrode})) || ...
    ~all(cellfun(@isscalar,{continuous_config(:).electrode})) || ...
    ~all(cellfun(@isreal,{continuous_config(:).electrode})) || ...
    ~all(cellfun(@isfinite,{continuous_config(:).electrode})) || ...
    ~all(cellfun(@(c) (c > 0),{continuous_config(:).electrode})) || ...
    ~all(cellfun(@(c) (round(c) == c),{continuous_config(:).electrode}))
  warning('electrode field must be a real positive integer scalar');
  bool = false;

elseif ~all(cellfun(@isnumeric,{continuous_config(:).channel})) || ...
    ~all(cellfun(@isscalar,{continuous_config(:).channel})) || ...
    ~all(cellfun(@isreal,{continuous_config(:).channel})) || ...
    ~all(cellfun(@isfinite,{continuous_config(:).channel})) || ...
    ~all(cellfun(@(c) (c >= 0),{continuous_config(:).channel})) || ...
    ~all(cellfun(@(c) (round(c) == c),{continuous_config(:).channel}))
  warning('channel field must be a non-negative integer scalar');
  bool = false;

elseif ~all(cellfun(@isnumeric,{continuous_config(:).depth})) || ...
    ~all(cellfun(@isscalar,{continuous_config(:).depth})) || ...
    ~all(cellfun(@isreal,{continuous_config(:).depth})) || ...
    ~all(cellfun(@isfinite,{continuous_config(:).depth})) || ...
    ~all(cellfun(@(c) (round(c) == c),{continuous_config(:).depth}))
  warning('depth field must be an integer scalar');
  bool = false;

elseif ~iscellstr({continuous_config(:).hemisphere}) || ...
    ~all(strcmp({continuous_config(:).hemisphere}, ...
    regexp({continuous_config(:).hemisphere},'[a-zA-Z0-9_]*','match','once')))
  warning(['hemisphere field must be a string containing characters ' ...
      '[a-zA-Z0-9_]']);
  bool = false;

elseif ~iscellstr({continuous_config(:).region}) || ...
    any(cellfun(@isempty,{continuous_config(:).region}))
  warning(['region field must be a non-empty string']);
  bool = false;

elseif ~all(cellfun(@isnumeric,{continuous_config(:).reference})) || ...
    ~all(cellfun(@isvector,{continuous_config(:).reference})) || ...
    ~all(cellfun(@(c) (size(c,1) == 1),{continuous_config(:).reference})) || ...
    ~all(cellfun(@(c) (size(c,2) == 2),{continuous_config(:).reference})) || ...
    ~all(cellfun(@(c) isreal(c),{continuous_config(:).reference})) || ...
    ~all(cellfun(@(c) (c(1) >= 0),{continuous_config(:).reference})) || ...
    ~all(cellfun(@(c) (c(2) >= 0),{continuous_config(:).reference}))
  warning(['reference field must be a 2-element row vector of real ' ...
      'non-negative integer values']);
  bool = false;

elseif ~all(cellfun(@isnumeric,{continuous_config(:).passband})) || ...
    ~all(cellfun(@isvector,{continuous_config(:).passband})) || ...
    ~all(cellfun(@(c) (size(c,1) == 1),{continuous_config(:).passband})) || ...
    ~all(cellfun(@(c) (size(c,2) == 2),{continuous_config(:).passband})) || ...
    ~all(cellfun(@(c) all(isreal(c)),{continuous_config(:).passband})) || ...
    ~all(cellfun(@(c) all(c > 0),{continuous_config(:).passband})) || ...
    ~all(cellfun(@(c) c(2) > c(1),{continuous_config(:).passband}))
  warning(['passband field must be a 2-element row vector of real ' ...
      'positive values in increasing order']);
  bool = false;

elseif ~all(cellfun(@isnumeric,{continuous_config(:).dspnum})) || ...
    ~all(cellfun(@isscalar,{continuous_config(:).dspnum})) || ...
    ~all(cellfun(@isreal,{continuous_config(:).dspnum})) || ...
    ~all(cellfun(@(c) ismember(c,DSPNUM_ALLOWED),{continuous_config(:).dspnum}))
  warning('dspnum field must be one of the allowed values');
  bool = false;

elseif ~all(cellfun(@isnumeric,{continuous_config(:).dspchan})) || ...
    ~all(cellfun(@isscalar,{continuous_config(:).dspchan})) || ...
    ~all(cellfun(@isreal,{continuous_config(:).dspchan})) || ...
    ~all(cellfun(@(c) ismember(c,DSPCHAN_ALLOWED),{continuous_config(:).dspchan}))
  warning('dspchan field must be one of the allowed values');
  bool = false;

elseif ~all(cellfun(@isnumeric,{continuous_config(:).Fs})) || ...
    ~all(cellfun(@(c) isa(c,'double'),{continuous_config(:).Fs})) || ...
    ~all(cellfun(@isscalar,{continuous_config(:).Fs})) || ...
    ~all(cellfun(@isreal,{continuous_config(:).Fs})) || ...
    ~all(cellfun(@isfinite,{continuous_config(:).Fs})) || ...
    ~all(cellfun(@(c) (c > 0),{continuous_config(:).Fs}))
  warning('Fs field must be a real positive scalar');
  bool = false;

% tricky here: nested cellfun calls!
elseif ~all(cellfun(@iscellstr,{continuous_config(:).sources})) || ...
    ~all(cellfun(@(c) numel(c) == 2,{continuous_config(:).sources})) || ...
    ~all(cellfun(@(c) (exist(c{1}) == 2) & (exist(c{2}) == 2), ...
    {continuous_config(:).sources})) || ...
    ~all(cellfun(@(c) isdir(fileparts(c{1})) & isdir(fileparts(c{2})), ...
    {continuous_config(:).sources}))
  warning(['sources field must be a 2-element cell array of strings, each ' ...
      'of which must specifies a filename, including full path']);
  bool = false;

else
  [i1, i2] = ndgrid(1:numel(continuous_config));
  pair_idx = find(i1 ~= i2);
  i1 = i1(pair_idx);
  i2 = i2(pair_idx);
  % Check that no two elements of continuous_config duplicate the same
  % (subject, day, electrode, channel) 4-tuple
  if any(arrayfun(@(s1,s2) isequal( ...
      {s1.subject, s1.day, s1.electrode, s1.channel}, ...
      {s2.subject, s2.day, s2.electrode, s1.channel}), ...
      continuous_config(i1),continuous_config(i2)))
    bool = false;
  % Check that no two elements of continuous_config which share the same 
  % (subject, day, electrode) 3-tuple have discrepant depth or region fields
  elseif any( arrayfun(@(s1,s2) isequal( ...
      {s1.subject, s1.day, s1.electrode}, ...
      {s2.subject, s2.day, s2.electrode}), ...
      continuous_config(i1),continuous_config(i2)) & ...
      arrayfun(@(s1,s2) ...
      ~isequal(s1.depth,s2.depth) || ~isequal(s1.region,s2.region), ...
      continuous_config(i1),continuous_config(i2)) )
    bool = false;
  else
    bool = true;
  end

end
