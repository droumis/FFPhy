% KERNLH: Likelihood of 2D kernal density function as a function of h, the 
%         global smoothing parameter, based on jackknifed density estimates.
%
%     Usage: rl = kernlh(h,xi)
%
%         h =     smoothing parameter (std of kernel function).
%         xi =    [n x 2] matrix of coordinates for observed points.
%         -------------------------------------------------------
%         rl =    reciprocal log-likelihood (to be minimized).
%

% Brunsdon, C. 1995. Estimating probability surfaces for geographic point data: 
%   an adaptive kernel algorithm. Computers & Geosciences 21:877-894.

% RE Strauss, 5/6/00

function rl = kernlh(h,xi)
  [n,p] = size(xi);

  fx = zeros(n,1);
  for i = 1:n                           % Jackknifed probability estimate
    xxi = xi;                           %   for each point with fixed h
    xxi(i,:) = [];
    fx(i) = kernden(xi(i,:),xxi,h);
  end;
  
  rl = -sum(log(fx));                   % Negative log-likelihood

  return;
