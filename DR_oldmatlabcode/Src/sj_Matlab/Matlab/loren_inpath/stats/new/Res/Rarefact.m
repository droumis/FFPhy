% RAREFACT: Calculation of expected number of species at a standard sample size, 
%           by the rarefaction method.  Randomized confidence intervals are 
%           estimated by randomizing the observed absolute frequencies; see 
%           prbcount().
%
%     Usage: [ES,ESci,S,N] = rarefact(freq,{iter},{CI_level},{doplot})
%
%           freq =      vector (length S) of absolute frequencies (counts) of 
%                         individuals per species, for S species.
%           iter =      optional number of randomization iterations at each value 
%                         of n, for confidence intervals [default = 0].
%           CI_level =  optional level of confidence level [default = 95].
%           doplot =    optional boolean flag indicating, if true, that 
%                         rarefaction curve is to be plotted [default = 0].
%           --------------------------------------------------------------------
%           ES =        [N x 1] vector of expected number of species in a random 
%                         sample of n individuals, for n = 1 to observed N.
%           ESci =      [N x 2] matrix of randomized confidence limits for 
%                         corresponding ES estimates.
%           S =         number of species in sample.
%           N =         total number of individuals in sample.
%

% Krebs, CJ. 1989. Ecological Methodology.  Harper & Row.  Chapter 10.
% Hurlbert, SH. 1971. The non-cept of species diversity: a critique and 
%   alternative parameters. Ecology 52:577-586.
% Sanders, H.L. 1968. Marine benthic diversity: a comparative study. Am. Nat. 
%   102:243-282.

% RE Strauss, 2/19/00

function [ES,ESci,S,N] = rarefact(freq,iter,CI_level,doplot)
  if (nargin < 2) iter = []; end;
  if (nargin < 3) CI_level = []; end;
  if (nargin < 4) doplot = []; end;

  esci_out = 0;
  if (nargout > 1)
    esci_out = 1;
  end;

  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(CI_level))
    CI_level = 95;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  if (~isvector(freq))
    error('  RAREFACT: frequencies must be in vector form.');
  end;

  freq = freq(:);
  i = find(freq==0);
  if (~isempty(i))
    freq(i) = [];
  end;

  S = length(freq);
  N = sum(freq);

  ES = zeros(N,1);
  ESci = [];

  for n = 1:N                             % Rarefaction curve
    es = 0;
    for s = 1:S
      es = es + (1-comb(N-freq(s),n)/comb(N,n));
    end;
    ES(n) = es;
  end;

  if (iter)                               % Randomized confidence intervals
    if (~esci_out & ~doplot)
      disp('  RAREFACT warning: no output, randomization not performed.');
      return;
    end;

    distrib = zeros(iter,N);                % Allocate sampling distribs
    specid = makegrps(1:S,freq);            % Vector of species identifiers

    for it = 1:iter                         % Iterate
      f = prbcount(freq./N,N,N,0,0);          % Randomize frequencies
      for n = 1:N                             % Rarefaction for each n
        es = 0;
        for s = 1:S
          es = es + (1-comb(N-f(s),n)/comb(N,n));
        end;
        distrib(it,n) = es;
      end;
    end;

    ESci = bootci(distrib,ES,[],CI_level)'; % Confidence intervals
  end;

  if (doplot)
    n = 1:N;
    plot(n,ES,'k');

    if (iter)
      hold on;
      plot(n,ESci(:,1),'k--');
      plot(n,ESci(:,2),'k--');
      putbnd([n,n,n]',[ES;ESci(:,1);ESci(:,2)]);
      hold off;
    else
      putbnd(n,ES);
    end;

    putxlab('Number of individuals in sample');
    putylab('Expected number of species');
  end;

  return;


