% CORRESP: Correspondence analysis of a contingency table of counts.
%          Returns subsets of row-eigenvectors and col-eigenvectors, along 
%          with percents variance explained.  
%          Optionally randomizes by generating random tables based on either
%          fixed or floating marginal totals.
%
%     Syntax: [row_crds,col_crds,percvar,CI_row_crds,CI_col_crds,CI_percvar] ...
%              = corresp(X,{nvect},{rand_iter},{marginals},{CI_level})
%
%          X =           [r x c] data matrix of counts.
%          nvect=        optional number of leading eigenvectors for which
%                          coordinates are desired [default=2].
%          rand_iter =   number of randomization iterations [default=0].
%          marginals =   TRUE if tables are to be randomized with fixed marginal
%                          totals, or FALSE with floating marginal totals
%                          [default=FALSE].
%          CI_level =    percent width of confidence intervals [default=95].
%          ---------------------------------------------------------------------
%          row_crds =    [r x nvect] matrix of row eigenvectors (columns) as
%                          coorespondence-axis coordinates.
%          col_crds =    [c x nvect] matrix of column eigenvectors (columns) as
%                          coorespondence-axis coordinates.
%          percvar =     [nvect x 1] vector of percents variance-explained
%                          for subsets of eigenvectors.
%          CI_row_crds = [r x 2*nvect] matrix of CI% asymmetric confidence
%                          limits of row coordinates, two columns per
%                          eigenvector (low,high).
%          CI_col_crds = [c x 2*nvect] matrix of CI% asymmetric confidence
%                          limits of row coordinates, two columns per
%                          eigenvector (low,high).
%          CI_percvar =  [nvect x 2] matrix of CI% asymmetric confidence limits
%                          of percents variance-explained.
%

% Everitt,BS & G Dunn. 1992. Applied Multivariate Data Analysis, p. 57-62.
%   Oxford Univ Press.

% RE Strauss, 11/28/95

function [row_crds,col_crds,percvar,CI_row_crds,CI_col_crds,CI_percvar] ...
            = corresp(X,nvect,rand_iter,marginals,CI_level)
  TRUE = 1; FALSE = 0;
  [r,c] = size(X);

  default_nvect = 2;                  % Defaults for input arguments
  default_rand_iter = 0;
  default_marginals = FALSE;
  default_CI_level = 95;

  if (nargin < 2)                     % If nvect null or not passed,
    nvect = [];                       %   set to default
  end;
  if (length(nvect) == 0)
    nvect = min([default_nvect,r,c]);
  end;

  if (nargin < 3)                     % Set missing randomization parameters
    rand_iter = default_rand_iter;    %   to defaults
  end;
  if (nargin < 4)
    marginals = default_marginals;
  end;
  if (nargin < 5)
    CI_level = default_CI_level;
  end;

  % Correspondence analysis
  total = sum(sum(X));                % Grand total of counts
  coltot = sum(X);                    % Marginal totals of probabilities
  rowtot = sum(X');

  exp = zeros(r,c);                   % Working matrix for expected values

  for i = 1:r                         % Sqrts of 'chi-square' cell contributions
    for j = 1:c
      exp(i,j) = rowtot(i)*coltot(j)/total;
      X(i,j) = (X(i,j)-exp(i,j))/sqrt(exp(i,j));
    end;
  end;

  [row_crds,sing_vals,col_crds] = svd(X); % Singular-value decomposition
  evals = diag(sing_vals .* sing_vals);
  percvar = 100 * evals / sum(evals); % Percents variance

  row_crds = row_crds(:,1:nvect);     % Retain subsets
  col_crds = col_crds(:,1:nvect);
  percvar = percvar(1:nvect);

  % Randomize correspondence analysis
  CI_row_crds = [];                  % Allocate return arguments
  CI_percvar = [];

  if (rand_iter > 0)
    clk = clock;                     % Seed random-number generator
    rand('seed',clk(6));             %   from system clock

    Bp = [percvar'];                 % Begin with initial solution
    for p=1:nvect                    % One Brow & Bcol matrix per eigenvector
         eval([ 'Brow' int2str(p) '=[row_crds(:,' int2str(p) ')''];' ]);
         eval([ 'Bcol' int2str(p) '=[col_crds(:,' int2str(p) ')''];' ]);
    end;
    rand_iter = rand_iter + 1;       % Include initial solution in count

    for bit = 1:rand_iter            % Bootstrap iterations
      Bx = continrn(rowtot,coltot,marginals); % Random contingency table
      Bx = (Bx-exp)./sqrt(exp);               % Sqrts of 'chi-square' cell contribs

      [Brow_crds,Bsing_vals,Bcol_crds] = svd(Bx); % Singular-value decomposition
      Bevals = diag(Bsing_vals .* Bsing_vals);
      Bpercvar = 100 * Bevals / sum(Bevals);  % Percents variance

      Brow_crds = Brow_crds(:,1:nvect);       % Retain subsets
      Bcol_crds = Bcol_crds(:,1:nvect);
      Bpercvar = Bpercvar(1:nvect);

      Bp = [Bp; Bpercvar'];                   % Accumulate results
      for p=1:nvect
        % Stash each eigenvector in a separate accum matrix
        eval([ 'Brow' int2str(p) '=[Brow' int2str(p) ...
               ';Brow_crds(:,' int2str(p) ')''];' ]);
        eval([ 'Bcol' int2str(p) '=[Bcol' int2str(p) ...
               ';Bcol_crds(:,' int2str(p) ')''];' ]);
      end;
    end;  % End bootstrap iterations

    high_limit = ceil((rand_iter*CI_level/100)); % CI indices
    low_limit = rand_iter - high_limit + 1;

    Bp = sort(Bp);                               % Sort variables independently
    CI_percvar =  Bp([low_limit,high_limit],:)'; % Confidence intervals

    for p=1:nvect
      eval([ 'Brow' int2str(p) '=sort(Brow' int2str(p) ');' ]);
      eval([ 'CI_row_crds(:,[2*' int2str(p) '-1,2*' int2str(p) ...
             '])=Brow' int2str(p) '([low_limit,high_limit],:)'';' ]);

      eval([ 'Bcol' int2str(p) '=sort(Bcol' int2str(p) ');' ]);
      eval([ 'CI_col_crds(:,[2*' int2str(p) '-1,2*' int2str(p) ...
             '])=Bcol' int2str(p) '([low_limit,high_limit],:)'';' ]);
    end;

    disp(sprintf('  Bootstrap time: %5.2f minutes\n', etime(clock,clk)/60));
  end;  % End bootstrap

  return;

