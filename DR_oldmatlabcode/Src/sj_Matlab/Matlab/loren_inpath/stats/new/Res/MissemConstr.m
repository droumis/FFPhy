% MISSEMCONSTR: Constrained maximum likelihood estimation of missing data, means, and
%         covariances for multivariate normal samples.  Estimated values are constrained
%         to within the interval [minest,maxest].  It is assumed that input data have 
%         already been log-transformed, if appropriate.
%
%     Usage: [X,M,C,propmiss,timeout,stop_time] = missemconstr(Y,{minest},{maxest},{maxtime})
%
%           Y =         [n x p] data matrix, including one or more missing 
%                         values indicated by values of NaN, Inf, or -Inf.
%           minest =    minimum allowable value for estimates; may be a scalar or a
%                         [1 x p] vector of values (one per variable).
%           maxest =    maximum allowable value for estimates; may be a scalar or a
%                         [1 x p] vector of values.
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

% RE Strauss, 1/7/02, modified from missem().


% Little, RJA & DB Rubin. 1987. Statistical Analysis with Missing Data. Wiley.

function [X,M,C,propmiss,timeout,stop_time] = missemconstr(Y,minest,maxest,maxtime)
  if (nargin < 2) minest = []; end;
  if (nargin < 3) maxest = []; end;
  if (nargin < 4) maxtime = []; end;

  tol = 1e-4;                         % Convergence tolerance for covariances

  [n,p] = size(Y);
  M = zeros(p,1);                     % Allocate return matrices
  C = zeros(p,p);
  timeout = 0;
  
  if (~isempty(minest))
    minest = minest(:)';
    if (isscalar(minest))
      minest = minest * ones(1,p);
    end;
  end;
  if (~isempty(maxest))
    maxest = maxest(:)';
    if (isscalar(maxest))
      maxest = maxest * ones(1,p);
    end;
  end;

  % Check for observations having all values missing

  for i = 1:n
    if (sum(finite(Y(i,:)))==0)
      error('  MISSEM: At least one observation has all values missing');
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
    C(i,i) = var(y);
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
        c = y'*y / (length(indx)-1);    % Adjusted covariance
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
                   
      if (~isempty(minest))               % Constrain estimates if requested
        if (any(any(X(i,varmiss) < ones(ni,1)*minest(varmiss))))
          for ivm = varmiss
            xx = X(i,ivm);
            ir = find(xx < minest(ivm));
            if (~isempty(ir))
              xx(ir) = minest(ivm)*ones(length(ir),1);
              X(i,ivm) = xx;
            end;
          end;
        end;
      end;
      if (~isempty(maxest))
        if (any(any(X(i,varmiss) > ones(ni,1)*maxest(varmiss))))
          for ivm = varmiss
            xx = X(i,ivm);
            ir = find(xx > maxest(ivm));
            if (~isempty(ir))
              xx(ir) = maxest(ivm)*ones(length(ir),1);
              X(i,ivm) = xx;
            end;
          end;
        end;
      end;
      
      c(varmiss,varmiss) = c(varmiss,varmiss) + F*ni;  % Adjust covars
    end; % for ip = 1:pattern

    % Little & Rubin, eqns 8.1, 8.2

    M = sum(X)'/n;                        % New estimates of means
    dev = X - ones(n,1)*M';               % Deviations from new means
    CP = dev'*dev;                        % Sufficient cross-products
    C = (CP+c)/(n-1);                     % New estimates of covariances

    Mdelta = max(abs(M-Msave));           % Check whether solution has converged
    Cdelta = max(max(abs(C-Csave)));

    if (max([Mdelta Cdelta]) < tol)
      do_again = 0;
    end;
  end; % while

  return;
