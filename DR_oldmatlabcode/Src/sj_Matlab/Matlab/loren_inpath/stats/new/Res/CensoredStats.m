% CensoredStats: Mean, standard deviation, their bootstrapped standard errors, and the
%                median for a singly truncated variable, truncated either above or below,
%                assuming that the variable is normally distributed.
%                Assumes that censored values are passed as NaN.  
%                    The mean and standard deviation are estimated using the method of 
%                "filling in with expected values", where the expected values are predicted 
%                from expected normal scores ("rankits").  Because this procedure strongly
%                assumes that the data are normally distributed, the data can be first 
%                transformed to near-normality using the Box-Cox transformation.  If the 
%                number of censored observations is less than the number of non-censored, 
%                then the median is within the range of the non-censored values; if not, 
%                it is extrapolated from predicted values.  Statistics are estimated by column.
%                    If no bootstrap iterations are indicated, then asymptotic estimates
%                of standard errors are provided based on numbers of non-censored observations.  
%                
%     Usage: [m,s,med,n,sem,ses,lambda,prob] = ...
%                                 censoredstats(X,{iter},{noplot},{noboxcox},{direction})
%
%         X =         [n x p] data matrix.
%         iter =      optional number of iterations for standard error estimate
%                       [default = 5000 if output arguments or the Box-Cox probabilities
%                       are requested].
%         noplot =    optional boolean flag indicating, if true, that a histogram of the
%                       observed and predicted data is not to be produced [default = 0].
%         noboxcox =  optional boolean flag indicating, if true, that the input data are not
%                       to be transformed to near-normality before estimating parameters
%                       [default = 0].
%         direction = optional direction of truncation:
%                       -1 (or other negative value) = truncated below [default];
%                       +1 (or other positive value) = truncated above.
%         ----------------------------------------------------------------------------------
%         m =         [1 x p] vector of means.
%         s =         corresponding vector of standard deviations.
%         med =       corresponding vector of medians.
%         n =         corresponding vector of non-censored sample sizes.
%         sem =       corresponding vector of standard errors of the means.
%         ses =       corresponding vector of standard errors of standard deviation.
%         lambda =    corresponding vector of Box-Cox parameters, if data are transformed.
%         prob =      corresponding vector of p-values for the null hypothesis that the 
%                       Box-Cox-transformed data are sampled from a normal distribution
%                       (based on 'iter' randomized iterations).  
%                       If this value is small (e.g., p<0.05), then the rankits procedure
%                       should not be used for estimating means and standard deviations. 
%

% RE Strauss, 10/16/02
%   10/21/02 - set default bootstrap iterations to 5000.
%   1/1/03 -   guard against bootstrapped samples with too few non-censored values;
%              revised sequence of output arguments; added output of lambda.
%   1/7/03 -   added normality test for Box-Cox-transformed data;
%              skip Box-Cox for data vectors having no NaNs.
%   2/9/03 -   don't do regression prediction when there are no missing data;
%              output lambda=NaN and prob=NaN if no missing data;
%              optional histgram of predicted and observed values.
%   6/11/03 -  change headsize in putarrow() to default value.

function [m,s,med,n,sem,ses,lambda,prob] = censoredstats(X,iter,noplot,noboxcox,direction)
  if (nargin < 1) help censoredstats; return; end;
  
  if (nargin < 2) iter = []; end;
  if (nargin < 3) noplot = []; end;
  if (nargin < 4) noboxcox = []; end;
  if (nargin < 5) direction = []; end;
  
  if (nargout < 5) iter = 0; end;
  if (nargout < 8)
    getprob = 0;
  else
    getprob = 1;
  end;
  
  if (isempty(iter))      iter = 5000; end;
  if (isempty(noplot))    noplot = 0; end;
  if (isempty(noboxcox))  noboxcox = 0; end;
  if (isempty(direction)) direction = -1; end;
  
  if (isvector(X))
    X = X(:);
  end;
  
  [n,p] = size(X);
  censored = (sum(~isfinite(X))>0);
  Xorig = X;
  
  lambda = NaN*ones(1,p);
  prob =   NaN*ones(1,p);
  Xp =     zeros(size(X));
  minX =   zeros(1,p);
  minXp =  zeros(1,p);
  
  if (~noboxcox)                          % Box-Cox transformation to normality
    c = zeros(1,p);
    for ip = 1:p
      x = X(:,ip);
      if (censored(ip))
        [x,lambda(ip),c(ip)] = boxcoxnorm(x,0,direction);
        X(:,ip) = x;
        if (getprob)
          [W,prob(ip)] = normaltest(x,iter);
        end;
      end;
    end;
  end;
  
  if (direction>0)                        % If data are right-censored, make left-sensored
    X = -X;
  end;

  for ip = 1:p
    i = find(isfinite(Xorig(:,ip)));
    minX = min(Xorig(i,ip));
    minXp = min(X(i,ip));
  end;
  
  m =   zeros(1,p);                       % Allocate output
  sem = zeros(1,p);
  s =   zeros(1,p);
  ses = zeros(1,p);
  med = zeros(1,p);
  nnc = zeros(1,p);
  
  for ip = 1:p                            % Cycle thru variables
    x = sort(X(:,ip));                      % Isolate and sort current var
    n = length(x);
    
    i = find(isfinite(x));                  % Count and delete censored values
    nmiss = n-length(i);
    nnc(ip) = nmiss;
    x = x(i);
    
    if (censored(ip))
      z = rankits(n);                       % Predict censored values from regression
      zx = z((nmiss+1):n);                  %   of data on normal scores (rankits)
      zpred = z(1:nmiss);
      [b,stats,pred] = linregr(zx,x,[],zpred);
      xp = [x;pred];
    else
      xp = x;
    end;
    Xp(:,ip) = xp;
    m(ip) = mean(xp);
    s(ip) = std(xp);
    med(ip) = median(xp);
    
    if (iter)                               % Boostrapped standard errors
      bootmean = zeros(iter,1);
      bootstd = zeros(iter,1);
      for it = 1:iter
        xb = 0;
        while (sum(isfinite(xb))<3)
          xb = sort(bootsamp(X(:,ip)));       % Bootstrap censored + non-censored values
        end;
        i = find(isfinite(xb));               % Count and delete censored values       
        nmiss = n-length(i);
        xb = xb(i);
        zx = z((nmiss+1):n);
        zpred = z(1:nmiss);
        [b,stats,pred] = linregr(zx,xb,[],zpred);
        xp = [xb;pred];
        bootmean(it) = mean(xp);
        bootstd(it) = std(xp);
      end;
      sem(ip) = std(bootmean);
      ses(ip) = std(bootstd);
    else
      sem(ip) = s(ip)/sqrt(n);
      ses(ip) = s(ip)/sqrt(2*n);
    end;
  end;
  
  if (direction>0)                          % If data are right-censored, adjust sign of mean
    m = -m;
    minX = -minX;
    minXp = -minXp;
  end;
  n = n - nnc;
 
  if (~noboxcox)                            % Invert Box-Cox transform and average the deviations from mean
    for ip = 1:p
      stats = [m(ip),m(ip)-sem(ip),m(ip)+sem(ip),...                                               % 1-3
               m(ip)-s(ip),m(ip)+s(ip),...                                                         % 4-5
               m(ip)-s(ip)-ses(ip),m(ip)-s(ip)+ses(ip),m(ip)+s(ip)-ses(ip),m(ip)+s(ip)+ses(ip),... % 6-9
               med(ip)];                                                                           % 10
      if (censored(ip))
        stats = real(boxcoxinv(stats,lambda(ip),c(ip)));
      end;
      m(ip) =   stats(1);
      sem(ip) = mean([abs(stats(1)-stats(2)),abs(stats(3)-stats(1))]);
      s(ip) =   mean([abs(stats(1)-stats(4)),abs(stats(5)-stats(1))]);;
      ses(ip) = mean([abs(stats(4)-stats(6)),abs(stats(7)-stats(4)),...
                      abs(stats(5)-stats(8)),abs(stats(9)-stats(5))]);
      med(ip) = stats(10);
    end;
  end;
  
%   if (~noplot & noboxcox)                 % Histogram of original data only
%     for ip = 1:p
%       xp = Xp(:,ip);
%       figure;
%       [nplot,freqs] = histgram(xp,[],[],[],[],1);
%       centers = freqs(:,2);
%       n = freqs(:,1);
%       shaftlen = 0.12;
%       v = axis;
%       [position,shaftlen,v] = histarw(m(ip),n,centers,v,shaftlen);
%       axis(v);
%       putarrow(position(3:4),position(1:2));
%       putxlab('Original data');
%       hold on;
%       delta = 0.004*(v(2)-v(1));
%       plot([minX(ip),minX(ip)],[0,0.97*v(4)],'k--');
%       plot([minX(ip)+delta,minX(ip)+delta],[0,0.95*v(4)],'w--');
%       hold off;
%       if (p>1)
%         puttitle(sprintf('Var %d',ip));
%       end;
%     end;
%   end;
  
  if (~noplot & ~noboxcox)                % Histogram of predicted and observed values
    for ip = 1:p
      xp = Xp(:,ip);
      
      figure;
      histgram(xp);
      putxlab('Transformed data');
      puttext(0.05,0.90,sprintf('lambda=%4.2f',lambda(ip)));
      v = axis;
      hold on;
      delta = 0.004*(v(2)-v(1));
      plot([minXp(ip),minXp(ip)],[0,0.97*v(4)],'k--');
      plot([minXp(ip)+delta,minXp(ip)+delta],[0,0.95*v(4)],'w--');
      hold off;
      if (p>1)
        puttitle(sprintf('Var %d - Box-Cox transformation',ip));
      else
        puttitle('Box-Cox transformation');
      end;
      
      xp = real(boxcoxinv(xp,lambda(ip),c(ip)));
      figure;
      [nplot,freqs] = histgram(xp,[],[],[],[],1);
      centers = freqs(:,2);
      nn = freqs(:,1);
      shaftlen = 0.12;
      v = axis;
      [position,shaftlen,v] = histarw(m(ip),nn,centers,v,shaftlen);
      axis(v);
%       putarrow(position(3:4),position(1:2),[],0.035);
      putarrow(position(3:4),position(1:2));
      putxlab('Original data');
      hold on;
      delta = 0.004*(v(2)-v(1));
      plot([minX(ip),minX(ip)],[0,0.97*v(4)],'k--');
      plot([minX(ip)+delta,minX(ip)+delta],[0,0.95*v(4)],'w--');
      hold off;
      if (p>1)
        puttitle(sprintf('Var %d',ip));
      end;
    end;
  end;

  return;
  