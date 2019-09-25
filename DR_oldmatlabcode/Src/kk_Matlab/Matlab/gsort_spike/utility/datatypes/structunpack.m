function names = structunpack(s, force)
%STRUCTUNPACK      Copies the fields of a structure into the workspace.
%   STRUCTUNPACK(S) creates a variable in the workspace with name
%   'fieldname' for each 'fieldname' that is a field of S, and copies the
%   value of S.('fieldname') into the variable 'fieldname'.
%
%   STRUCTUNPACK(S,1) will overwrite existing variables of name
%   'fieldname' without warning.  STRUCTUNPACK(S,0) is equivalent to
%   STRUCTUNPACK(S) and will issue a warning that 'fieldname' exists but
%   will not overwrite the existing variable.
%
%   FIELDS = STRUCTUNPACK(S, ...) returns the fieldnames that were
%   successfully copied to the workspace.

if (nargin < 2), force = 0; end;

names = fieldnames(s);
for n = 1:length(names)
	if (~force && (exist(names{n}) == 1))  % if not overwriting, detect name conflicts
		warning(['Field ' names{n} ' already exists in workspace.  Will not copy.']);
		names(n) = '';
	else
		assignin('caller', names{n}, s.(names{n}));
	end
end

if (nargout == 0), clear names; end;
