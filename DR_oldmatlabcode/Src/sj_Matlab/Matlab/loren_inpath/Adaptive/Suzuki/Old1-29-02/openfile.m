function f = openfile (filename, mode)

% openfile.m
%
% like fopen, but it aborts execution on failure
%
% Last modified: 26 Nov 98

[f, msg] = fopen (filename, mode);
if f == -1
   error (msg);
end 