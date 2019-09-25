function ds = pushd(dir)
% pushd - change directory and remember old location

global DIRSTACK
dirs;

if (nargin == 0)
  if (length(DIRSTACK) == 0)
    error ('no previous directory on stack');
  end
  
  thisd = pwd;
  popd;
  DIRSTACK = {thisd, DIRSTACK{:}};
end

if (nargin == 1)
  DIRSTACK = {pwd, DIRSTACK{:}};
  cd(dir);
end

ds = dirs;
