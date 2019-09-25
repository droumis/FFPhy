% HAUSDORF: Estimates the fractal Hausdorff dimension for a path
%
%     Usage: [D,se,r2,Lmax,pathlen] = hausdorf(path,{nsteps},{noplot})
%
%           path =    [n x 2] matrix of point coordinates for a 2D path.
%           step =    optional vector of step numbers used to estimate Hausdorff
%                       dimension [default = n:2n].
%           noplot =  optional boolean vector indicating that plots of the path and 
%                       regression are to be suppressed [default = 0].
%           ----------------------------------------------------------------
%           D =       Hausdorff dimension.
%           se =      standard error of D.
%           r2 =      regression r-squared.
%           Lmax =    asymptotic length of path.
%           pathlen = exact length of path.
%

% Estimate of Lmax is incorrect.
% Must introduce nonlinear regression to find the fractal dimension, and consider 
%   Reeve (1992) about the standard error.
% See Slice (1992) for discussion of how to determine which step-lengths to use.  Modify 
%   for use of step-lengths directly (log-distributed) rather than step numbers?

function [D,se,r2,Lmax,pathlen] = hausdorf(path,nsteps,noplot)
  [n,p] = size(path);

  if (nargin < 2)                               % Check for input arguments
    nsteps = [];
  end;
  if (nargin < 3)
    noplot = [];
  end;

  if (isempty(nsteps))                          % Default values
    nsteps = [n:2*n]';
  end;
  if (isempty(noplot))
    noplot = 0;
  end;

  if (size(nsteps,1)==1)                        % Transform nsteps to col vector
    nsteps = nsteps';
  end;

  if (nargout > 2)
    c1 = path(1:(n-1),:);
    c2 = path(2:n,:);
    len = sqrt( (c1(:,1)-c2(:,1)).^2 + ...
                (c1(:,2)-c2(:,2)).^2 );
    pathlen = sum(len);
  end;

  steplen = zeros(size(nsteps));
  totlen =  zeros(size(nsteps));

  for s = 1:length(nsteps)
    steplen(s) = pathseg(path,nsteps(s));
    totlen(s) =  steplen(s) * nsteps(s);
  end;

  log_steplen = log(steplen);
  log_totlen = log(totlen);

  [b,stats] = linregr(log_steplen,log_totlen);
  s2 = stats(2);
s2
  x = log_steplen - mean(log_steplen);
  se = sqrt(s2./(x'*x));
se
  se = se * 1.5*(length(nsteps)-2);           % Reeve (1992) correction

  D = 1-b(2);
  Lmax = exp(b(1)./D);                        % <-- this is wrong

  if (~noplot)
    figure;
    plot(path(:,1),path(:,2),'o');
    hold on;
    plot(path(:,1),path(:,2));
    hold off;
    putbnd(path(:,1),path(:,2));

    x = [min(log_steplen) max(log_steplen)];
    y = b(2)*x + b(1);

    figure;
    plot(log_steplen,log_totlen,'o');
    hold on;
    plot(x,y);
    hold off;
    putbnd(log_steplen,log_totlen);
    putxlab('Ln(Step length)');
    putylab('Ln(Path length)');

    x = linspace(steplen(1),steplen(length(steplen)));
    y = exp(b(2)*log(x) + b(1));

    figure;
    plot(steplen,totlen,'o');
    hold on;
    plot(x,y);
    hold off;
    putbnd(steplen,totlen);
    putxlab('Step length');
    putylab('Path length');
  end;

  return;
