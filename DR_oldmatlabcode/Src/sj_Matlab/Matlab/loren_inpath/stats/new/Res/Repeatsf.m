% REPEATSF: Objective function for repeatsf().
%
%     Usage: [totx2,expfreq] = repeatsf(monofreq,seq,obsfreq,monoseq)
%
%           monofreq = [4 x 1] vector of estimated mono-nucleotide frequencies.
%           seq =      [r x maxlen] matrix of nucleotide-sequence identifiers 
%                        (numeric equivalents of ascii identifiers).
%           obsfreq =  corresponding observed absolute frequencies (counts).
%           monoseq =  [4 x 1] vector of unique identifiers (numeric).
%           -------------------------------------------------------------------
%           totx2 =    total chi-square deviation of expected from observed.
%           expfreq =  vector of expected frequencies.
%

% RE Strauss, 2/19/98
%   9/20/99 - update handling of null input arguments.

function [totx2,expfreq] = repeatsf(monofreq,seq,obsfreq,monoseq)
  [r,maxseqlen] = size(seq);

  expfreq = ones(size(obsfreq));

  for nuc = 1:4                         % Expected relative frequencies
    [i,j] = find(seq == monoseq(nuc));
    if (~isempty(i))
      expfreq(i) = expfreq(i)*monofreq(nuc);
    end;
  end;
  expfreq = expfreq * sum(obsfreq);     % Expected absolute frequencies

  % Objective function

  dev = obsfreq - expfreq;
  expfreq = max([expfreq'; 1e-6*ones(size(expfreq))'])';
  totx2 = sum(dev.*dev./expfreq);       % Chi-square deviations

  % Constraints

%  constr(1) = abs(sum(monofreq)-1);     % Freqs must sum to unity
%  constr(2) = -min(monofreq-0.01);      % All freqs must be >0
%  constr(3) = -totx2;                   % Objective function must be >0

  return;



