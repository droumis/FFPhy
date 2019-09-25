function bool = is_session(session)
%IS_SESSION Validate whether the input is an SESSION struct array
%
%   IS_SESSION(SESSION) returns true if SESSION is a struct array that
%   contains valid sessions meta-data.
%
%   See also READ_SESSION, written by smk.
%
%Depends on:
%   IS_TIMESTRING (written by smk)
%   IS_TIMERANGE (written by smk)
%   
% Written by SMK 2009 June 22.
%

if (exist('is_timestring') ~= 2)
  error('IS_SESSION depends on m-file IS_TIMESTRING (written by smk)');
end
if (exist('is_timerange') ~= 2)
  error('IS_SESSION depends on m-file IS_TIMERANGE (written by smk)');
end

REQUIRED_FIELDS = { ...
    'subject'     , ...
    'day'         , ...
    'epoch'       , ...
    'environment' , ...
    'tstart'      , ...
    'tend'        , ...
    'timerange'   , ...
    'source'      };

if isempty(session) && isstruct(session) && ...
    all(isfield(session,REQUIRED_FIELDS))
  warning('SESSION struct is empty');
  bool = true;

elseif ~isstruct(session) || ~all(isfield(session,REQUIRED_FIELDS))
  warning('missing required field(s)');
  bool = false;

elseif ~iscellstr({session(:).subject}) || ...
    ~all(strcmp({session(:).subject}, ...
    regexp({session(:).subject},'[a-zA-Z0-9]*','match','once')))
  warning(['subject field must be a string containing characters ' ...
      '[a-zA-Z0-9]']);
  bool = false;

elseif ~all(cellfun(@isnumeric,{session(:).day})) || ...
    ~all(cellfun(@isscalar,{session(:).day})) || ...
    ~all(cellfun(@isreal,{session(:).day})) || ...
    ~all(cellfun(@isfinite,{session(:).day})) || ...
    ~all(cellfun(@(c) (c >= 0),{session(:).day})) || ...
    ~all(cellfun(@(c) (round(c) == c),{session(:).day}))
  warning('day field must be a real non-negative integer scalar');
  bool = false;

elseif ~iscellstr({session(:).epoch}) || ...
    ~all(strcmp({session(:).epoch}, ...
    regexp({session(:).epoch},'[a-zA-Z0-9]*','match','once')))
  warning(['epoch field must be a string containing characters ' ...
      '[a-zA-Z0-9]']);
  bool = false;

elseif ~iscellstr({session(:).environment}) || ...
    ~all(strcmp({session(:).environment}, ...
    regexp({session(:).environment},'[a-zA-Z0-9]*','match','once')))
  warning(['environment field must be a string containing characters ' ...
      '[a-zA-Z0-9]']);
  bool = false;

elseif ~all(cellfun(@isnumeric,{session(:).timerange})) || ...
    ~all(cellfun(@(c) isa(c,'uint32'),{session(:).timerange})) || ...
    ~all(cellfun(@isvector,{session(:).timerange})) || ...
    ~all(cellfun(@(c) (size(c,2) == 2),{session(:).timerange})) || ...
    ~all(cellfun(@(c) (c(2) > c(1)),{session(:).timerange}))
  warning(['timerange field must be a 2-element row vector of uint32 ' ...
      'timestamps in increasing order']);
  bool = false;

elseif ~all(cellfun(@is_timestring,{session(:).tstart}))
  warning('tstart field must be a valid time string');
  bool = false;

elseif ~all(cellfun(@is_timestring,{session(:).tend}))
  warning('tend field must be a valid time string');
  bool = false;

elseif ~all(arrayfun(@(s)s.timerange(1)==str2ts(s.tstart),session(:))) || ...
    ~all(arrayfun(@(s)s.timerange(2)==str2ts(s.tend),session(:)))
  warning('tstart and tend fields must match timerange field');
  bool = false;

elseif ~iscellstr({session(:).source}) || ...
    ~all(cellfun(@(c) (exist(c) == 2),{session(:).source})) || ...
    ~all(cellfun(@(c) isdir(fileparts(c)), {session(:).source}))
  warning(['source field must be a string that specifies a fileepoch, ' ...
      'including path']);
  bool = false;

else

  % this is a sentinel value which will be overwritten if the remaining tests
  % fail
  bool = true;

  for subject = unique({session(:).subject});
    per_subject = session(strcmp(subject,{session(:).subject}));
    for day = unique([per_subject.day])
      per_subject_per_day = per_subject([per_subject.day] == day);
      if numel(unique({per_subject_per_day(:).epoch})) ~= ...
          numel(per_subject_per_day)
        warning(['session data structs for the same subject and day ' ...
            'have duplicate epoch designations']);
        bool = false;
      end
      if ~is_timerange(sortrows(vertcat(per_subject_per_day(:).timerange)))
        warning(['session data structs for the same subject and day ' ...
            'have overlapping time ranges']);
        bool = false;
      end
    end
  end

end

