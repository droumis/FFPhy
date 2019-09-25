function tf = struct_cmp(s1,s2,varargin)
%STRUCT_CMP Compare selected fields of two structs.
%
%   TF = STRUCT_CMP(S1,S2,FIELDS), where S1 is a scalar struct and S2 can be
%   a scalar struct or a struct array, and FIELDS is a cell array of strings,
%   performs the following comparison:
%   
%   cellfun(@(f) isequal(s1.(FIELDS{1}),f),{s2.(FIELDS{1})) &
%   cellfun(@(f) isequal(s1.(FIELDS{2}),f),{s2.(FIELDS{2})) &
%   ...
%   cellfun(@(f) isequal(s1.(FIELDS{end}),f),{s2.(FIELDS{end}))
%
%   for each element in FIELDS. An error if raised if S1 or S2 is missing one of
%   the fieldnames specified in FIELDS. The return value TF is a logical array
%   of the same size as S2.
%
%   STRUCT_CMP(S1,S2) is equivalent to STRUCT_CMP(S1,S2,fieldnames(S1))
%
%Written by SMK 2009 November 9
%

if (nargin > 3)
  error('too many arguments');
end

if (nargin == 3)
  fields = varargin{1};
elseif (nargin == 2)
  if ~isstruct(s1) || ~isscalar(s1)
    error('S1 must be a scalar struct');
  end
  fields = fieldnames(s1);
end

if ~iscellstr(fields) || isempty(fields)
  error('fields must be a non-empty cell array of strings');
end
if ~isstruct(s1) || ~all(isfield(s1,fields))
  error('S1 is not a struct with the required fields');
end
if ~isscalar(s1)
  error('S1 must be a scalar struct');
end
if ~isstruct(s2) || ~all(isfield(s2,fields))
  error('S2 is not a struct with the required fields');
end

if isempty(s2)
  tf = false(size(s2));
  return;
end

tf = true([1 numel(s2)]);
for i = 1:numel(fields)
  tf = tf & cellfun(@(f) isequal(s1.(fields{i}),f),{s2.(fields{i})});
end
tf = reshape(tf,size(s2));


