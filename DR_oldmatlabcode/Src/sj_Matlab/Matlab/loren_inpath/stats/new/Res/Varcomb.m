% VARCOMB: Given an [n x p] data matrix with missing data, finds the number of 
%           complete observations for all possible combinations of 1,...,p 
%           characters.  Practical for up to 20 or so variables.  Returns the 
%           ln(condition) or ln(reciprocal condition).  The condition factors are
%           defined in function condfactor().
%
%     Usage: [nobs,nmiss,combvars,combobs,condfact,toomany,C] = ...
%              varcomb(X,{condtype},{minobs},{minvars},{usecorr},{sortcond},{printwarn})
%
%           X =         [n x p] data matrix.
%           condtype =  optional type of condition factor to be used:
%                         1 = condition factor, which measures stability of the 
%                             eigensolution: log((eigenvalue 1)/(eigenvalue 2));
%                         0 = difference between condition factor of submatrix and
%                             original matrix [default];
%                        -1 = reciprocal condition factor, which measures 
%                             stability under matrix inversion:
%                             log((eigenvalue 2)/(eigenvalue 1)).
%           minobs =    optional minimum number of observations to be considered 
%                         [default = 3].
%           minvars =   optional minimum number of variables in combination 
%                         [default = 3].
%           usecorr =   optional boolean variable indicating, if true, that the 
%                         conditions of the correlation matrices rather than 
%                         covariance matrices are to be evaluated [default = 0].
%           sortcond =  optional boolean flag indicating, if true, that the 
%                         combinations of observations and variables are to be 
%                         sorted by decreasing condition [default = 1].
%           printwarn = optional boolean flag indicating, if true, that none-fatal 
%                         warning messages are to be printed for error conditions
%                         [default = 0].
%           ---------------------------------------------------------------------
%           nobs =      vector of numbers of observations.
%           nmiss =     corresponding matrix of numbers of missing values in 
%                         complete submatrix.
%           combvars =  corresponding matrix of combinations of variables.  
%                         For each value of 'nobs', only the largest subsets 
%                         of variables are provided.
%           combobs =   corresponding matrix of combinations of observations.
%           condfact =  condition factors for corresponding covariance matrices 
%                         (used to sort output, high to low).
%           toomany =   boolean flag indicating that there are too many missing values.
%           C =         cov/corr matrix of original data matrix, adjusted to be
%                         positive-definite if necessary.
%

% RE Strauss, 3/9/99
%   9/29/00 -  if more than a threshold number of vars, reduce to subset having 
%                most observations.
%   11/2/00 -  added option to sort by decreasing condition factor.
%   11/11/00 - added choice of condition factor;
%              changed maximum number of variables from 15 to 19;
%              changed minimum number variables from 2 to 3.
%   11/12/00 - changed 'sortcond' default to 1.
%   11/13/00 - return condition factors as ln(condition) rather than condition.
%   12/18/00 - changed 'minobs' default from p+1 to 3.
%   5/4/01 -   added list of included observations to output.
%   7/16/01 -  added 'matchorig' option.
%   7/17/01 -  change condition factors to use first two eigenvalues (in 
%                decreasing magnitude) rather than first and last.
%   5/14/02 -  added failure flag for too many missing values.
%   5/17/02 -  add flag for optional warning messages;
%              changed finite() to isfinite().
%   5/30/02 -  implement error flag for call to posdef().
%   6/21/02 -  return 'C' as output argument.
%   7/1/02 -   isolated some code as covpairwise().


function [nobs,nmiss,combvars,combobs,condfact,toomany,C] = ...
           varcomb(X,condtype,minobs,minvars,usecorr,sortcond,printwarn)
  if (nargin<2) condtype = []; end;
  if (nargin<3) minobs = []; end;
  if (nargin<4) minvars = []; end;
  if (nargin<5) usecorr = []; end;
  if (nargin<6) sortcond = []; end;
  if (nargin<7) printwarn = []; end;

  maxvar = 19;                          % Max permissible number of vars

  if (isempty(minobs))          % Default input arguments
    minobs = 3;
  end;
  if (isempty(minvars))
    minvars = 3;
  end;
  if (isempty(condtype))
    condtype = 0;
  end;
  if (isempty(usecorr))
    usecorr = 0;
  end;
  if (isempty(sortcond))
    sortcond = 1;
  end;
  if (isempty(printwarn))
    printwarn = 0;
  end;
  
  switch (condtype)             % Expand condtype argument
    case -1,
      matchorig = 0;
      condtype = 0;
    case 0,
      matchorig = 1;
      condtype = 1;
    case 1,
      matchorig = 0;
      condtype = 1;
    otherwise
      error('  VARCOMB: invalid condition type.');
  end;
  
  toomany = 0;

  [n,p] = size(X);

  if (p>=maxvar)
    if (printwarn)
      disp(['  VARCOMB warning: number of variables reduced to ',int2str(maxvar),...
          ' having the']);
      disp( '                   largest number of observations.');
    end;

    b = isfinite(X);
    sumb = sum(b);
    [sumb,i] = sort(-sumb);
    i = sort(i(1:maxvar));
    X = X(:,i);
    if (printwarn)
      disp('  Variables kept:');
      disp(i);
    end;

    [n,p] = size(X);
  end;

  if (matchorig)                        % If want to match original matrix,
    [C,M,N] = covpairwise(X,1);         %   estimate covariance structure of original matrix
    
    if (any(any(N==0)))
      if (printwarn)
        if (any(diag(N)==0))
          disp('  VARCOMB warning: too many missing values for 1 or more variables');
        else
          disp('  VARCOMB warning: 1 or more pairs of variables produce no covariances');
        end;
      end;
      nobs = [];
      nmiss = [];
      combvars = [];
      combobs = [];
      condfact = [];
      toomany = 1;
      return;
    end;

    condfactorig = condfactor(C,condtype);
  end;

  XX = isfinite(X);               % Transform data matrix to boolean

  if (sum(sum(~XX)) == 0)         % If no missing values, quit
    nobs = n;
    nmiss = 0;
    combvars = [1:p];
    combobs = [1:n];
    condfact = NaN;
    return;
  end;

  nc = sum(combin(p,minvars:p));
  nobs_full =  zeros(nc,1);
  nmiss_full = zeros(nc,1);
  combvars_full = zeros(nc,p);

  i = 0;
  for r = p:-1:minvars          % Cycle backwards through sets of vars
    c = combvals(p,r);            % List of combinations of variables
    for ic = 1:size(c,1)          % Cycle thru list
      cc = c(ic,:);
      x = XX(:,cc);                   % Extract vars
      ntot = sum(prod(x')');          % Number of obs for this combination
      nmiss = sum(~x(:));             % Number of missing values for this combination
      i = i+1;
      nobs_full(i) = ntot;
      nmiss_full(i) = nmiss;
      combvars_full(i,1:length(cc)) = cc;
    end;
  end;

  nobs = [];                        % Restrict output to largest sets of vars
  combvars = [];                       %   for each unique value of nobs
  combobs = [];
  nmiss = [];

  nvar = sum((combvars_full>0)')';

  for ni = n:-1:minobs          
    i = find(nobs_full == ni);
    if (~isempty(i))
      nobs_red = nobs_full(i);
      nmiss_red = nmiss_full(i);
      nvar_red = nvar(i);
      ncombvars_red = combvars_full(i,:);
      maxnvar = max(nvar_red);
      i = find(nvar_red==maxnvar);
      nobs =  [nobs; nobs_red(i)];
      combvars = [combvars; ncombvars_red(i,:)];
      nmiss = [nmiss; nmiss_red(i)];
    end;
  end;

  nvars = sums(combvars);           % Compact combination matrix to get rid of cols of zeros
  i = find(nvars>0);
  combvars = combvars(:,i);
  nc = size(combvars,1);

  condfact = zeros(nc,1);
  for ci = 1:nc                     % Get covar-matrix condition factors for combinations
    i = find(combvars(ci,:)>0);        %   of vars and observations
    cc = combvars(ci,i);
    x = X(:,cc);                      % Extract vars from data matrix
    s = sum(x')';                     % Row sums
    i = find(isfinite(s));            % Obs w/o missing values
    x = x(i,:);
    combobs(ci,1:length(i)) = i';     % Save list of obs for this combination of vars

    if (usecorr)                      
      cmat = corrcoef(x);
    else
      cmat = cov(x);
    end;
    
    condfact(ci) = condfactor(cmat,condtype);
  end;

  if (matchorig)                    % Optionally get deviation from original cond factor
    condfact = abs(condfactorig - condfact);
  end;

  if (sortcond)                     % Sort matrices
    if (matchorig)                    % By increasing deviation from orig cond factor
       [condfact,nobs,nmiss,combvars,combobs] = ...
                                  sortmat(condfact,nobs,nmiss,combvars,combobs);
    else                              % Or by decreasing condition factor
      [condfact,nobs,nmiss,combvars,combobs] = ...
                                  sortmat(-condfact,nobs,nmiss,combvars,combobs); 
      condfact = -condfact;
    end;
  end;
  
  if (isempty(nobs) | isempty(nmiss) | isempty(combvars) | isempty(combobs) | isempty(condfact))
    toomany = 1;
  end;

  return;

