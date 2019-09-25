% MISSEM: Maximum likelihood estimation of missing data, means, and
%         covariances for multivariate normal samples.  
%         It is assumed that input data have already been log-transformed, 
%         if appropriate.
%
%     Usage: [X,M,C,propmiss,timeout,stop_time] = missem(Y,{maxtime})
%
%           Y =         [n x p] data matrix, including one or more missing 
%                         values indicated by values of NaN, Inf, or -Inf.
%           maxtime =   optional number of seconds of elapsed time to allow 
%                         before terminating with an error message [default = 
%                         no maximum].
%           -------------------------------------------------------------------
%           X =         input data matrix with imputed missing values.
%           M =         [p x 1] vector of means.
%           C =         [p x p] covariance matrix.
%           propmiss =  proportion of missing data in Y.
%           timeout =   boolean flag indicating that a timeout (total time 
%                         exceeding 'maxtime') has occurred; if so, other output 
%                         matrices are set to null.
%           stop_time = clock time at which timeout occurred.
%

% RE Strauss, 3/14/96
%   7/9/99 -  convergence criterion changed.
%   4/10/00 - return proportion of missing data.
%   4/17/00 - added terminal error messages.
%   10/2/00 - added maximum-time option.
%   1/25/01 - return stop time if timeout occurs.
%   4/12/02 - find complex modulus of any complex numbers.
%   4/10/03 - use N rather than N-1 to estimate covariances.

% Little, RJA & DB Rubin. 1987. Statistical Analysis with Missing Data. Wiley.

function [X,M,C,propmiss,timeout,stop_time] = missem(Y,maxtime)
  if (nargin < 2) maxtime = []; end;

  tol = 1e-4;                         % Convergence tolerance for covariances

  [n,p] = size(Y);
  M = zeros(p,1);                     % Allocate return matrices
  C = zeros(p,p);
  timeout = 0;

  % Check for observations having all values missing

  for i = 1:n
    if (sum(finite(Y(i,:)))==0)
      error('  MISSEM: at least one observation has all values missing');
    end;
  end;

  % Initial values of parameter set: estimate means and variances from all
  % available values of each variable, and estimate covariances from all
  % available pairwise values, using means as above (Little & Rubin, eqn 3.5).

  for i = 1:p                         % For all variables,
    y = Y(finite(Y(:,i)),i);            % Eliminate missing values
    if (isempty(y))
      error('  MISSEM: too many missing values for 1 or more variables');
    end;
    M(i) =  mean(y);                    % Estimate means & variances
%     C(i,i) = var(y);
    leny = length(y);
    C(i,i) = var(y)*(leny-1)/leny;
  end;

  nc = 0;
  sumc = 0;

  for i = 1:(p-1)                     % For all pairs of variables,
    for j = (i+1):p
      yi = finite(Y(:,i));              % Find pairwise available values
      yj = finite(Y(:,j));
      indx = find(yi & yj);
      y = Y(indx,[i j]);
      y(:,1) = y(:,1) - M(i);           % Subtract respective means
      y(:,2) = y(:,2) - M(j);
      if (length(indx)>1)
%         c = y'*y / (length(indx)-1);    % Adjusted covariance
        c = y'*y / length(indx);    % Adjusted covariance
        sumc = sumc + c;
        nc = nc + 1;
        C(i,j) = c(1,2);
        C(j,i) = c(1,2);
      else
        error('  MISSEM: 1 or more pairs of variables produce no covariances');
      end;
    end;
  end;

  [i,j] = find(C==0);                   % If any covariances still zero,
  if (length(i)>0)                      %   set to mean covariance
    meanc = sumc / nc;
    for k = 1:length(i)
      C(i(k),j(k)) = meanc;
    end;
  end;

  % Find observations having identical patterns of missing values

  misspatn = zeros(n,1);
  missing = ~finite(Y);
  obsmiss = find(~finite(sum(Y')));
  misspatn(obsmiss) = 10e6 * ones(size(obsmiss));

  propmiss = sum(sum(missing))/(n*p); % Proportion of missing data

  pattern = 0;
  for i=1:n                           % Number the patterns
    if (misspatn(i)>pattern)            % Get next observation
      pattern = pattern+1;
      misspatn(i) = pattern;
      if (i<n)                          % If not last observation,
        for j = (i+1):n                 %   search for other obs with same pattern
          if (misspatn(j)>pattern & all(missing(i,:)==missing(j,:)))
            misspatn(j) = pattern;
          end;
        end;
      end;
    end;
  end;

  % Converge to solution

  t0 = clock;
  stop_time = [];

  do_again = 1;
  while (do_again)
    if (~isempty(maxtime))
      if (etime(clock,t0) > maxtime)
        stop_time = clock
        disp('  MISSEM: time limit exceeded.');
        X = [];
        M = [];
        C = [];
        propmiss = [];
        timeout = 1;
        return;
      end;
    end;

    Msave = M;                          % Save current parameter set
    Csave = C;
    X = Y;                              % Copy input matrix, with missing values
    c = zeros(p,p);                     % Matrix of var/covar adjustments

    for ip = 1:pattern                  % Est missing vals for obs for each pattern
      i = find(misspatn==ip);
      ni = length(i);
      varhave = find(finite(Y(i(1),:)));  % List of observed vars for this obs
      nhave = length(varhave);
      varmiss = find(~finite(Y(i(1),:))); % List of missing vars for this obs
      [D,E,F] = sweepreg(M,C,varhave);    % Regress missing vars on observed vars
      X(i,varmiss) = ones(ni,1)*E(1,:) ...             % Estimate missing values
                     + Y(i,varhave)*E(2:(nhave+1),:);  %   and plug into X
      c(varmiss,varmiss) = c(varmiss,varmiss) + F*ni;  % Adjust covars
    end;

    % Little & Rubin, eqns 8.1, 8.2

    M = sum(X)'/n;                        % New estimates of means
    dev = X - ones(n,1)*M';               % Deviations from new means
    CP = dev'*dev;                        % Sufficient cross-products
%     C = (CP+c)/(n-1);                     % New estimates of covariances
    C = (CP+c)/n;                     % New estimates of covariances

%    if (any(Msave<tol))                   % Check whether solution has converged
      Mdelta = max(abs(M-Msave));    
%    else
%      Mdelta = max(abs(M-Msave)./Msave);    
%    end;

%    if (any(Csave<tol))
      Cdelta = max(max(abs(C-Csave)));
%    else
%      Cdelta = max(max(abs(C-Csave)./Csave));
%    end;

%    delta = [Mdelta Cdelta]
    if (max([Mdelta Cdelta]) < tol)
      do_again = 0;
    end;
  end; % while
  
  if (any(imag(X(:))))
    [i,j] = find(imag(X));
    for k = 1:length(i)
      X(i(k),j(k)) = abs(X(i(k),j(k)));
    end;
  end;

  return;
