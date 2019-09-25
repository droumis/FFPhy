% KERNDEN: Kernel density functions f(x) for 2D x, given a set of observed 2D 
%          data points xi and corresponding values of h, the smoothing 
%          parameter.  The kernal density estimator is bivariate normal.
%
%     Usage: fx = kernden(x,xi,h)
%
%         x =   [m x 2] matrix of coordinates for points at which function 
%                 is to be evaluated.
%         xi =  [n x 2] matrix of observed point coordinates.
%         h =   smoothing parameter (std of kernel function), either a scalar 
%                 (global value) or a vector of values corresponding to the xi.
%         ---------------------------------------------------------------------
%         fx =  [m x 1] vector of density estimates corresponding to x.
%

% Worton, BJ. 1989. Kernel methods for estimating the utilization distribution 
%   in home-range studies. Ecology 70:164-168.
% Brunsdon, C. 1995. Estimating probability surfaces for geographic point data: 
%   an adaptive kernel algorithm. Computers & Geosciences 21:877-894.

% RE Strauss, 5/7/00

function fx = kernden(x,xi,h)
  [m,pm] = size(x);
  [n,pn] = size(xi);

  fx = zeros(m,1);                        % Allocate output matrix

  for i = 1:m                             % Evaluate at each point
    d = (ones(n,1)*x(i,:)-xi)';
    d = sqrt(sum(d.*d))';                   % Dists of x(i) from the observed pts
    f = normpdf(d,0,h);                     % Prob at each distance
    fx(i) = mean(f);                        % Mean of all probs
  end;

  return;
