function dirs = pushd(dir)
% popd - change back to old directory

global DIRSTACK

if (nargin == 0)
  if (length(DIRSTACK) == 0)
    error ('no previous directory on stack');
  end

  cd (DIRSTACK{1});
  DIRSTACK = {DIRSTACK{2:end}};
end

dirs = char(DIRSTACK);
