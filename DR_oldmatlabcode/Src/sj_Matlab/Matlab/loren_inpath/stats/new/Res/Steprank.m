% STEPRANK: To cure a singular correlation or covariance matrix with rank r < p 
%           for p variables, uses a stepwise procedure to find the subset of 
%           r variables such that the log-ratio of the first two eigenvalues of the 
%           reduced matrix is similar as possible to that of the original matrix.
%
%     Usage: vars = steprank(C,{r})
%
%           C =    square symmetric matrix.
%           r =    optional number of reduced variables [default = rank(C)].
%           -------------------------------------------------------------------
%           vars = row vector of indices of the subset of variables (unsorted).
%

% RE Strauss, 11/25/98
%   9/27/01 - use issqsym() to check matrix structure;
%             compare condition  of reduced matrix against that of full matrix.
%   7/2/02 -  call condfactor() for condition factor;
%             allow input of r.

function vars = steprank(C,r)
  if (nargin < 2) r = []; end;

  [n,p] = size(C);

  if (~issqsym(C))
    error('  STEPRANK: input matrix must be square symmetric.');
  end;

  if (isempty(r))
    r = rank(C);                      % Rank of original matrix
    if (r < 2)
      error('  STEPRANK: input matrix is of rank 1');
    end;
  else
    if (r < 2)
      error('  STEPRANK: at least two variables must be returned.');
    end;
  end; 
  
  v = [1:p];
  if (r == p)
    vars = v;
    return;
  end;

  orig_cond = condfactor(C);          % Condition of original matrix
  if (~isfinite(orig_cond))
    vars = 1:r;
    disp('  STEPRANK warning: covariance matrix ill-conditioned; variable subset not optimal.');
    return;
  end;

  best_diff = Inf;
  for i = 1:(p-1)                     % Find two best vars
    for j = (i+1):p
      cf = condfactor(C([i,j],[i,j]));  % Condition of reduced matrix
      diff = abs(cf-orig_cond);
      if (diff < best_diff)
        best_diff = diff;
        isave = i;
        jsave = j;
      end;
    end;
  end;

  vars = [isave jsave];               % Initialize lists of variables
  v(vars) = [];

  if (r > 2)                          % Find vars 3 ... r
    for ir = 3:r
      best_diff = Inf;
      for i = 1:length(v)
        j = [vars v(i)];
        bc = condfactor(C(j,j));        % Condition of reduced matrix
        diff = abs(bc-orig_cond);
        if (diff < best_diff)
          best_diff = diff;
          isave = i;
        end;
      end;
      vars = [vars v(isave)];
      v(isave) = [];
    end;
  end;

  return;

