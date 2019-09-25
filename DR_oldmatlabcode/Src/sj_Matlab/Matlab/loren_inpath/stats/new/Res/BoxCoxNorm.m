% BoxCoxNorm: A modified Box-Cox transformation to find the best power 
%         transformation of an empirical distribution to normality by maximizing 
%         the Shapiro-Francia W' statistic, which is approximately the squared 
%         correlation between the sorted data and the corresponding rankits or 
%         normal scores.  Assumes that missing data are left- or right-censored data 
%         (passed as NaN's). If any of the data are non-positive, a constant is added 
%         to the vector to make them positive.
%
%     Syntax:  [xp,lambda,c,W,Wp] = boxcoxnorm(x,{plotflag},{censdir})
%
%         x =        vector of observations for a single variable.
%         plotflag = optional boolean flag indicating, if true, that a histogram of
%                      the transformed data and superimposed normal distribution are 
%                      to be plotted [default = 0].
%         censdir =  optional boolean variable indicating the direction of censoring:
%                      -1 = left-censored [default];
%                      +1 = right-censored.
%         ----------------------------------------------------------------------------
%         xp =       corresponding vector of transformed variable.
%         lambda =   Box-Cox parameter.
%         c =        value added to data before transforming to ensure all positive 
%                      values.
%         W =        Shapiro-Francia W' statistic for original data.
%         Wp =       Shapiro-Francia W' statistic for transformed data.
%

% RE Strauss, 1/1/03, modified from boxcox().
%   1/3/03 - remove scaling of transformed data; 
%            return 'c'; 
%            add error for non-vector input.

function [xp,lambda,c,W,Wp] = boxcoxnorm(x,plotflag,censdir)
  if (nargin < 1) help boxcoxnorm; return; end;
  
  if (nargin < 2) plotflag = []; end;
  if (nargin < 3) censdir = []; end;
  
  if (isempty(plotflag)) plotflag = 0; end;
  if (isempty(censdir))  censdir = -1; end;
  
  censored = 0;
  if (any(~isfinite(x)))
    censored = 1;
  end;
  
  if (~isvector(x))
    error('  BoxCoxNorm: input data must be in vector form');
  end;
  
  i = find(isfinite(x));
  xmin = min(x(i));
  xmax = max(x(i));
  if (xmin <= 0)                        % Passed values must be positive
    c = abs(xmin)+1;
    x = x+c;
  else
    c = 0;
  end;
  
  xmean = mean(x(i));
  xstd = std(x(i));

  lambda = fminbnd('boxcoxnormf',-10,10,[],x,censdir);   % Optimize lambda
  
  xp = x;
  if (abs(lambda) > eps)                
    xp(i) = ((x(i).^lambda)-1)/lambda;
  else
    xp(i) = log(x(i));
  end;

  [W,prob,a,xa] = normaltest(x,0,censdir);      % Test-statistic value for original data
  [Wp,prob,ap,xap] = normaltest(xp,0,censdir);  % Test-statistic value for transformed data

  if (plotflag)
    scatter(xa,a);
    putregrline(xa,a);
    putxlab('Data');
    putylab('Coefficients');
    puttitle('Originial data');
    
    scatter(xap,ap);
    putregrline(xap,ap);
    putxlab('Data');
    putylab('Coefficients');
    puttitle('Transformed data');
    
    figure;
    histgram(x(i)-c);
    putxlab('Data');
    puttitle('Original data');
    
    figure;
    histgram(xp(i)-c);
    putxlab('Data');
    puttitle('Transformed data');
  end;

  return;

