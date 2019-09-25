% KWCRIT: Randomized critical values for the Kruskal-Wallace test, given the 
%         sample sizes.  Critical values are right-tailed.
%
%     Usage: [critval,X2val,alpha] = kwcrit(n,{nties},{alpha},{iter},{doplot})
%
%         n =       vector (length k) of sample sizes for k groups.
%         nties =   number of tied data values [default = 0].
%         alpha =   vector of alpha levels for which critical values are desired 
%                   [default = 0.05; 0.01; 0.005; 0.001].
%         iter =    number of randomization iterations [default = 1/min(alpha)].
%         doplot =  boolean flag indicating that histogram of null distribution 
%                     is to be produced [default = 0 = false].
%         ----------------------------------------------------------------------
%         critval = vector (of same size as 'alpha') of critical values.
%         X2val =   corresponding critical values from the chi-square 
%                     distribution with k-1 df.
%         alpha =   corresponding alpha levels.
%

% RE Strauss, 10/27/99

function [critval,X2val,alpha] = kwcrit(n,nties,alpha,iter,doplot)
  if (nargin < 2) nties = []; end;
  if (nargin < 3) alpha = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) doplot = []; end;

  if (isempty(nties))
    nties = 0;
  end;
  if (isempty(alpha))
    alpha = [0.05; 0.01; 0.005; 0.001];
  end;
  min_iter = 1/min(alpha);
  if (isempty(iter))
    iter = min_iter;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  if (iter < min_iter)
    iter = min_iter;
    disp('  KWCRIT warning: iterations too few to estimate critical value');
    disp('                  for smallest alpha.  Iterations reset to');
    disp(iter);
  end;

  critval = zeros(size(alpha));           % Allocate output matrix
  H = zeros(iter,1);                      % Allocate test-statistic distribution

  k = length(n);                          % Number of groups
  N = sum(n);                             % Total sample size
  g = makegrps(1:k,n);                    % Grouping vector
  r = 1:N;                                % Ranks

  for it = 1:iter                         % Produce randomized distribution
    if (nties>0)                            % Introduce ties, if specified
      r = 1:N;
      tiegrp = zeros(1,N);

      if (nties==1)
        tiepos = ceil(rand*(N-1));            % Random tie position
        rtie = r(tiepos)+0.5;
        r(tiepos:tiepos+1) = [rtie rtie];
      else
        rp = randperm(N-1);                   % Random tie positions
        tiepos = rp(1:nties);     
        for i = 1:nties                       % Scatter ties
          tp = tiepos(i);
          if (all(tiegrp(tp:tp+1)==[0 0]))
            tiegrp(tp:tp+1) = [i i];
          elseif (tiegrp(tp)>0)
            tiegrp(tp+1) = tiegrp(tp);
          else
            tiegrp(tp) = tiegrp(tp+1);
          end;
        end;

        for i = 1:nties                       % Find mean ranks
          j = find(tiegrp==i);
          if (~isempty(j))
            r(j) = mean(r(j))*ones(1,length(j));
          end;
        end;
      end;
    end;

    r = r(randperm(N));                     % Randomize ranks
    H(it) = kwstat(r,g,1);
  end;

  H = sort(H);
  pos = (1-alpha)*iter;
  for i = 1:length(alpha)
    if (isintegr(pos(i)))
      critval(i) = H(pos(i));
    else
      fp = floor(pos(i));
      cp = fp+1;
      critval(i) = H(fp) + (H(cp)-H(fp))*(pos(i)-fp);
    end;
  end;

  X2val = chi2inv(1-alpha,k-1);

  if (doplot)
    figure;
    histgram(H);
  end;

  return;
