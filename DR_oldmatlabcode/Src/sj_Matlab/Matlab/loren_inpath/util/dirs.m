function ds = dirs
% dirs - list directory stack

global DIRSTACK
if (isempty(DIRSTACK))
  DIRSTACK = {};
end

ds = char({pwd, DIRSTACK{:}});
