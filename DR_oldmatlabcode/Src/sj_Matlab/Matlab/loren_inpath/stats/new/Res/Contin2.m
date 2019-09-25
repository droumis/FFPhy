% CONTIN2:  Two-way contingency-table analysis, asymptotic or randomized.
%           If randomized, probabilities are based either on fixed
%           marginal totals (type 1) or floating marginal totals (type 2).
%
%     Syntax: [pr,totx2,df,cell_x2,cell_prob,exp] =
%                 contin2(obs,{useG},{iter},{fixed},{doplot})
%
%           obs =       [r x c] matrix of observed counts.
%           useG =      optional flag indicating whether or not the G-statistic 
%                         is to be used in place of the conventional X^2 
%                         [default = 0].
%           iter =      optional number of randomization iterations [default = 0].
%           fixed =     optional boolean flag indicating whether marginal (row & 
%                         column) totals to be treated as fixed [=1, the default] 
%                         or unconstrained [=0] for randomization.
%           doplot =    optional boolean flag specifying plot of overall randomized
%                         chi-square distribution [default = 0].
%           -----------------------------------------------------------------------
%           pr =        overall significance-level of observed table.
%           totx2 =     overall observed chi-square value.
%           df =        degrees of freedom for table.
%           cell_x2 =   [r x c] matrix of observed chi-square values, by cell.
%           cell_pr =   [r x c] matrix of chi-square probabilities, by cell.
%           exp =       [r x c] matrix of expected values.
%

% Romesberg,HC & K Marshall. 1985. CHITEST: a Monte-Carlo computer program for
%   contingency table tests.  Computers & Geosciences 11:69-78.

% RE Strauss, 11/11/97
%   5/14/99 - miscellaneous improvements.
%   2/28/00 - make histogram relative rather than absolute.
%   6/6/00 -  check for zero row/col totals.

function [pr,totx2,df,cell_x2,cell_pr,exp] = contin2(obs,useG,iter,fixed,doplot)
  if (nargin<2) useG = []; end;
  if (nargin<3) iter = []; end;
  if (nargin<4) fixed = []; end;
  if (nargin<5) doplot = []; end;

  if (isempty(useG))                      % Default input arguments
    useG = 0;
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(fixed))
    fixed = 1;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  if (~fixed & ~iter)
    error('  CONTIN2: unconstrained solution possible only by randomization.');
  end;

  [r,c] = size(obs);
  if (min([r,c])==1)
    error('  CONTIN2: use 1-way analysis (CONTIN1) for vector of observations.');
  end;

  exp =  zeros(r,c);                      % Allocate working matrices

  coltot = sum(obs);                      % Marginal totals
  rowtot = sum(obs')';
  N = sum(sum(obs));                      % Grand total

  if (any(coltot==0) | any(rowtot==0))
    error('  CONTIN2: at least one row or column total is zero.');
  end;

  exp = rowtot*coltot/N;                  % Expected values

  if (useG)                               % Observed cell G or chi-square values
    cell_x2 = -2.*exp.*log(obs./exp);       % Max-likelihood statistic
  else
    cell_x2 = (obs-exp).^2./exp;            % Pearson chi-square statistic
  end;

  totx2 = sum(sum(cell_x2));              % Observed total chi-square
  if (useG)                               % William's correction for overall G
    q = 1 + ((r*c+1))/(6*N);
    totx2 = totx2 ./ q;
  end;

  df = (r-1)*(c-1);                       % Degrees of freedom

  if (~iter)                              % Asymptotic chi-square probabilities
    cell_pr = 1-chi2cdf(cell_x2,df);        % By cell
    pr = 1-chi2cdf(totx2,df);               % Global

  else                                    % Randomized probabilities
    pr = 0;
    cell_pr = zeros(r,c);
    incr = 1/iter;
    if (doplot)
      randx2 = zeros(iter,1);
    end;

    for it = 1:iter                         % Iterate random table construction
      obsit = continrn(rowtot,coltot,fixed);  % Random table

      if (useG)                               % Observed cell G or chi-square values
        cell_x2it = -2.*exp.*log(obsit./exp);   % Max-likelihood statistic
      else
        cell_x2it = (obsit-exp).^2./exp;        % Pearson chi-square statistic
      end;

      totx2it = sum(sum(cell_x2it));            % Observed total chi-square
      if (useG)                               % William's correction for overall G
        totx2it = totx2it ./ q;
      end;

      if (doplot)                           % Optionally save for plot of distribution
        randx2(it) = totx2it;
      end;

      [i,j] = find(cell_x2it >= cell_x2);     % Cell probabilities
      if (~isempty(i))
        for ic = 1:length(i)
          cell_pr(i(ic),j(ic)) = cell_pr(i(ic),j(ic))+incr;
        end;
      end;

      if (totx2it >= totx2)                   % Total probability
        pr = pr + incr;
      end;
    end;  % Iterations

    if (doplot)                               % Plot randomized distribution
      histgram(randx2,[],[],[],[],[],'rel');
      xlabel('Randomized chi-square statistic');
    end;
  end;  % Randomization

  return;
