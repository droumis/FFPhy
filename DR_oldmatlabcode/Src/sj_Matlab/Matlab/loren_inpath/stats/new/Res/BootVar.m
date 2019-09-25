% BootVar: create a bootstrapped sampling distribution of the variance.
%
%     Usage: [var_x,stderr_b,stderr_t,CI] = bootvar(x,{iter},{CI_level})
%
%         x = [n x 1] vector of data.
%         iter = optional number of bootstrap iterations [default = 1000].
%         CI_level = optional level of confidence interval [default = 95].
%         ---------------------------------------------------------------------
%         var_x = variance of x.
%         stderr_b = standard error of variance from the sampling distribution.
%         stderr_t = theoretical standare error of the variance.
%         CI = 2-element vector containing lower and upper confidence limits.
%

function [var_x,stderr_b,stderr_t,CI] = bootvar(x,iter,CI_level)
  if (nargin < 2)
    iter = [];
  end;
  if (nargin < 3)
    CI_level = [];
  end;
  
  if (isempty(iter))
    iter = 1000;
  end;
  if (isempty(CI_level))
    CI_level = 95;
  end;
  
  x = x(:);                         % Force dataa to column vector
  var_x = var(x);                   % Variance of x
  
  distrib = zeros(iter,1);          % Allocate sampling distribution
  
  for it = 1:iter                   % Create the sampling distribution
    bx = bootsample(x);               % Bootstrapped sample
    distrib(it) = var(bx);            % Stash the variance
  end;
  
  histgram(distrib);                % Histogram of the distribution
  putxlab('Variance of X');
  
  stderr_t = var(x)*sqrt(2/length(x));  % Theoretical standard error
  stderr_b = std(distrib);            % Bootstrapped standard error
  
  L = (100-CI_level)/2;
  U = 100-L;
  CI = prctile(distrib,[L,U]);
 
  return;
  