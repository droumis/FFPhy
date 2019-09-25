% PROBIT: Maximum likelihood probit regression of binary response on continuous 
%         'dose' variable.
%           If the {0,1} responses do not overlap on the independent variable, 
%         there is no unique solution; a solution is obtained by assuming the 
%         gap to represent the central 99% of the normal distribution of response 
%         (tolerance distribution).
%
%     Usage: [b,ndp,maxmin,loglike] = probit(x,y,{doplot})
%
%         x =       [n x 1] vector of independent variable.
%         y =       [n x 1] vector of binary dependent variable.
%         doplot =  optional boolean flag indicating that plot is to be produced.
%         ------------------------------------------------------------------------
%         b =       vector of estimates of coefficients: [b1 b0].
%         ndp =     vector of mean and standard deviation of normal tolerance 
%                     distribution; the mean is the x-value of the inflection 
%                     point, and the std is a measure of variability of response.
%         maxmin =  [1 x 2] vector giving the maximum x at which no response was 
%                     observed, and the minimum x at which a response was observed.
%         loglike = log-likelihood of solution.
%

% RE Strauss, 4/2/99
%   9/3/99 -  transpose input row vectors to col vectors;
%             changed plot colors for Matlab v5.
%   1/4/00 -  changed fminu() to fmins().
%   6/19/01 - introduce gap if response/no-response ranges abutt.

function [b,ndp,maxmin,loglike] = probit(x,y,doplot)
  if (nargin<3) doplot=[]; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  [n,px] = size(x);
  [ny,py] = size(y);

  if (n==1 & px>1)
    x = x';
    [n,px] = size(x);
  end;
  if (ny==1 & py>1)
    y = y';
    [ny,py] = size(y);
  end;

  if (ny~=n | px>1 |py>1)
    error('  PROBIT: input vectors not compatible');
  end;

  if (min(y)~=0 | max(y)~=1 | length(unique(y))~=2)
    error('  PROBIT: dependent variable must be binary with values {0,1}');
  end;

  xmin = min(x);
  xmax = max(x);

  max0 = max(x(y==0));
  min1 = min(x(y==1));
  maxmin = [max0 min1];

  x = [x ones(n,1)];            % Append vector of ones to independent variables

  if (max0 == min1)             % If response/no-response abutt, create a gap
    max0 = max0-eps;
    min1 = min1+eps;
  end;
    
  if (max0 < min1)              % If response/no-response don't overlap, finish
    ndp = [(max0+min1)/2, (min1-max0)/2*norminv(0.995)];
    b = [1/ndp(2) -ndp(1)/ndp(2)];
    loglike = -probitf(b,x,y);

  else                           % Else if overlap,
    ndp0 = [(max0+min1)/2, (max0-min1)/2*norminv(0.975)];
    b0 = [1/ndp0(2) -ndp0(1)/ndp0(2)];  % Initial parameter estimates

    b = fmins('probitf',b0,[],[],x,y);  % Final max-likelihood estimates
    loglike = -probitf(b,x,y);          % Log-likelihood

    ndp = [-b(2)/b(1) abs(1/b(1))];     % Params of tolerance distribution
  end;

  if (doplot)
    x = x(:,1);
    xpred = linspace(xmin,xmax)';
    ypred = normcdf(b(1)*xpred+b(2));
    ypdf =  normpdf(b(1)*xpred+b(2));

    figure;
    plot(x,y,'ok');
    putbnd(x,y);
    hold on;
    plot(xpred,ypred,'k');
%    plot(xpred,ypdf,'b');
    hold off;
  end;

  return;
