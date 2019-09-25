% DISCRDET: Find the determinant of the between-to-within matrix from 
%           a discriminant analysis, as a function of the number of variables.
%
%     Syntax: [d,rnk,nv] = discrdet(X,grps,numvar,{repl})
%
%         X =      [n x p] data matrix (obs x vars).
%         grps =   row or column vector of group identifiers.
%         numvar = vector of numbers of random variables to include in 
%                    discriminant analyses.
%         repl =   number of randomized replicate analyses for each number of 
%                    variables [default = 100].
%         ---------------------------------------------------------------------------
%         d =      vector of mean log10-determinants corresponding to numbers 
%                    of variables.
%         rnk =    vector of mean matrix ranks corresponding to numbers of variables.
%         nv =     vector of numbers of variables (same as numvar).
%

% RE Strauss, 11/26/96
%   11/29/99 - changed calling sequence.

function [d,rnk,nv] = discrdet(X,grps,numvar,repl)
  if (nargin < 4) repl = []; end;

  if (isempty(repl))
    repl = 100;
  end;

  d = numvar;
  rnk = numvar;
  nv = numvar;

  index = uniquef(grps);
  ngrps = length(index);          % Number of groups
  [nobs,nvars] = size(X);         % Numbers of observations & variables
  G = design(grps);               % ANOVA-type design matrix
  vincl = randperm(nvars);        % Random variable list

  for n = 1:length(numvar)
    detrep = zeros(repl,1);       % Allocate vectors for replicated results
    rnkrep = zeros(repl,1);

    for rp = 1:repl               % Randomization replicates
      vincl = vincl(randperm(nvars)); % List of vars to include
      Xr = X(:,vincl(1:nv(n)));     % Reduced data matrix

      mean_W = (G'*G)\G'*Xr;        % Within-group means
      devs_W = Xr - G*mean_W;       % Within-group deviations
      W = devs_W'*devs_W;           % Within-group SSCP matrix

      devs_T = Xr - ones(nobs,1)*mean(Xr);  % Total deviations
      T = devs_T'*devs_T;           % Total SSCP matrix

      B = T - W;                    % Among-sample SSCP matrix
      detrep(rp) = det(B);
%      rnkrep(rp) = rank(B);
      rnkrep(rp) = rcond(B);
    end;
    d(n) = mean(abs(log10(detrep))) * sign(mean(detrep)-1);
    rnk(n) = mean(rnkrep);
  end;

  return;
