% WALKRAND: Models a random walk in p dimensions from the origin of the 
%           coordinate system.
%
%     Usage: [dist,walk] = walkrand({p},{nsteps},{stepflag})
%
%           p =        number of dimensions [default = 2].
%           nsteps =   number of steps      [default = 100].
%           stepflag = flag indicating distribution from which step-lengths 
%                       are sampled:
%                         0 = normal,  abs(N(0,1)) [default]
%                         1 = uniform, U[0,1]
%                         2 = constant 1
%           ---------------------------------------------------------------
%           dist =     Euclidean distance from origin at end of walk.
%           walk =     [nsteps+1 x p] matrix of path coordinates.
% 

% RE Strauss, 8/12/95

function [dist,walk] = walkrand(p,nsteps,stepflag)
  if (nargin < 1) p = []; end;
  if (nargin < 2) nsteps = []; end;
  if (nargin < 3) stepflag = []; end;

  const_step = 1;
  walk = zeros(nsteps+1,p);

  if (isempty(p))                     % Defaults
    p = 2;
  end;
  if (isempty(nsteps))
    nsteps = 100;
  end;
  if (isempty(stepflag))
    stepflag = 0;
  end;

  for step = 2:(nsteps+1)
    if (step==0)                      % Assign step length
      steplen = abs(randn);             % Random normal
    elseif (step==1)
      steplen = rand;                   % Random uniform
    elseif (step==2)
      steplen = const_step;             % Constant
    end;

    v = rand(1,p)*2-1;                % Assign random directional vector
    vlen = sqrt(v*v');                  % Length of resultant vector
    prop = (steplen/vlen)^2;            % Length-adjustment factor
    v = sign(v) .* sqrt(prop*v.^2);     % Length-adjusted vector

    walk(step,:) = walk(step-1,1) + v;  % Take step in p dimensions
  end;

  dist = eucl(walk(1,:),walk(nsteps+1,:)); % Dist from origin to last point

  return;

