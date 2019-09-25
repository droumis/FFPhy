% DIVERSITY: Species-diversity indices, estimated from a vector of relative or 
%            absolute frequencies.  Zero frequencies are removed before diversity 
%            is estimated.  Confidence intervals are estimated by treating the 
%            observed proportions (relative frequencies) and total individual 
%            count (N) as constants, and randomly varying the allocation of 
%            individuals to species; see function prbcount().
%
%         Diversity indices:
%         1) Shannon-Wiener index (H', using natural logs) = uncertainty of 
%            predicting the species of the next individual collected, given 
%            the current sample; range log(N/(N-S)) to log(S).
%         2) exp(H') (MacArthur) = number of equally common species that would 
%            produce the same diversity as H'; range N/(N-S) to S.
%
%         3) Simpson's index (D) = probability of picking two organisms at random 
%            that are of the same species (actually a measure of homogeneity 
%            rather than diversity); range 1/S to 1. Uses Pielou's (1969) finite-
%            population correction for absolute frequencies.
%         4) Complement of Simpson's index (1-D) = probability of picking two 
%            organisms at random that are of different species; range 0 to 1-1/S.
%         5) Reciprocal of Simpson's index (1/D) (Williams) = number of equally 
%            common species required to generate the observed heterogeneity of 
%            the sample; range 1 to S for relative frequencies, or 
%            1 to S(S+1)/2 for absolute frequencies.
%
%         Evenness indices:
%           For Shannon-Wiener indices: E = H'/(log S).
%           For Simpson indices:
%             Relative frequencies: E = (1-D)/(1-1/S).
%             Absolute frequencies: E = (1-D)/(1-Dmax), where Dmax is the 
%               diversity of the most evenly-spaced distribution possible of 
%               individuals among species.
%           
%
%     Usage: [divers,evenness,S,dCI,eCI] = ...
%                                 diversity(freqs,{kind},{N},{iter},{CI_level})
%
%         freqs =     [S x C] matrix of relative or absolute frequencies for S 
%                       species and C assemblages.
%         kind =      integer 1-5 (as described above) indicating the particular 
%                       diversity measure to be returned [default = 1].
%         N =         optional vector (length C) of total numbers of individuals 
%                       (required for the randomization of relative frequencies).
%                       If N is a scalar but C>1, N is assumed to be identical 
%                       for all semblages.
%         iter =      optional number of randomization iterations for confidence 
%                       intervals [default = 0].
%         CI_level =  optional level for confidence intervals [default = 95].
%         -----------------------------------------------------------------------
%         divers =    [1 x C] vector of diversity measures.
%         evenness =  [1 x C] vector of evenness measures.
%         S =         [1 x C] vector of numbers of species.
%         dCI =       [2 x C] matrix of lower and upper confidence limits for 
%                       diversity indices.
%         eCI =       [2 x C] matrix of lower and upper confidence limits for 
%                       evenness indices.
%

% Krebs, CJ. 1989. Ecological Methodology.  Harper & Row.  Chapter 10.

% RE Strauss, 2/17/00

function [divers,evenness,S,dCI,eCI] = diversity(freqs,kind,N,iter,CI_level)
  if (nargin < 2) kind = []; end;
  if (nargin < 3) N = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) CI_level = []; end;

  if (isempty(kind))                      % Default input arguments
    kind = 1;
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(CI_level))
    CI_level = 0.95;
  end;

  get_evenness = 0;
  if (nargout > 1)
    get_evenness = 1;
  end;

  if (isin(kind,[1 2]))
    kindname = 'shannon';
  elseif (isin(kind,[3,4,5]))
    kindname = 'simpson';
  else
    error('  DIVERSITY: invalid kind of diversity index specified');
  end;

  use_counts = 0;
  if (isintegr(freqs))                    % If absolute counts
    N = sum(freqs);                       %   Find N
    if (kindname == 'simpson')            %   Set flag for Simpson's index
      use_counts = 1;
    end;
  end;

  if (isvector(freqs))                    % Convert vector of freqs to col vector
    freqs = freqs(:);
  end;
  C = size(freqs,2);                      % Number of assemblages

  if (isempty(N))
    N = zeros(1,C);
  end;
  if (length(N)==1 & C>1)                 % Convert scalar N to vector
    N = N*ones(1,C);
  end;

  [divers,evenness,S] = diversf(freqs,kind,kindname,N,use_counts); % Get estimates

  dCI = [];                     
  eCI = [];

  if (iter)                               % Randomized confidence intervals
    if (isempty(N))
      error('  DIVERSITY: need sample sizes (N) to randomize relative frequencies');
    end;

    dsamp = zeros(iter,C);                % Allocate sampling distributions
    esamp = zeros(iter,C);
    f = freqs;

    for it = 1:iter                       % Randomized sampling distributions
      for ic = 1:C                              % For each assemblage,
        ff = freqs(:,ic);                       %   isolate frequencies
        ff = ff./N(ic);                         %   convert to proportions
        f(:,ic) = prbcount(ff,N(ic),N(ic),1,0); %   randomize counts, given proportions
      end;
      [dsamp(it,:),esamp(it,:)] = diversf(f,kind,kindname,N,use_counts); % Get estimates
    end;

    dCI = bootci(sort(dsamp),divers,[],CI_level);     % Confidence limits
    eCI = bootci(sort(esamp),evenness,[],CI_level);
  end;

  return;


