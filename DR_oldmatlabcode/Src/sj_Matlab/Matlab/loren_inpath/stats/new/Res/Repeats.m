% REPEATS:  Test the statistical significance of observed frequencies of tandom repeats 
%           (di-, tri-, tetra-nucleotides, etc.) by estimating the mononucleotide 
%           frequencies.  Finds the mono frequencies that minimize the difference 
%           between observed and expected frequencies, based on the conventional 
%           'chi-square' statistic.
%
%     Usage: [monoseq,monofreq,totx2,totprob,x2,prob,df] = repeats(seq,obsfreq)
%
%           seq =      [r x maxlen] matrix of nucleotide-sequence identifiers 
%                        (left-justified, in signel quotes, padded with blanks 
%                        on the right).
%           obsfreq =  corresponding observed absolute frequencies (counts); the 
%                        complementary nucleotides are assumed to be present in the same 
%                        frequencies.
%           ------------------------------------------------------------------------
%           monoseq =  [4 x 1] list of nucleotide identifiers.
%           monofreq = [4 x 1] corresponding estimated frequencies.
%           totprob =  overall significance-level of observed table.
%           totx2 =    overall observed value of chi-square statistic.
%           x2 =       [r x c] matrix of observed chi-square contributions, by cell.
%           prob =     [r x c] matrix of randomized chi-square probabilities, by cell
%           df =       degrees-of-freedom for test
%
%
%           Example of matrices of sequence identifiers and frequencies:
%
%               seq = [ 'GT  '        obsfreq = [ 190
%                       'CT  '                    110
%                       'TA  '                     76
%                       'GTT '                     14
%                       'GAT '                     18
%                       'GTA '                     26
%                       'GATA'                      8
%                       'GTAT'                      4
%                       'GGAA' ]                   11 ]
%

% RE Strauss, 2/19/98
 
function [monoseq,monofreq,totx2,totprob,x2,prob,df] = repeats(seq,obsfreq)
  [r,maxseqlen] = size(seq);

  monoseq = 'ACGT';             % Nucleotide bases
  monoseq = abs(monoseq);       % Convert to numeric values

  seq = abs(seq);               % Convert ascii identifiers to numbers

  f = linspace(0.01,0.49,100);  % Vary f(G) freq from 0.01 - 0.49
  ff = f;                       % F(f) = total chi-square, to be minimized

  for i = 1:length(f)           % Cycle thru possible frequencies
    fG = f(i);
    fC = fG;
    fT = 0.5 - fG;
    fA = fT;

    monofreq = [fA fC fG fT];
    ff(i) = repeatsf(monofreq,seq,obsfreq,monoseq);
  end;

  [totx2,minpos] = min(ff);
  fG = f(minpos);
  fC = fG;
  fT = 0.5 - fG;
  fA = fT;
  monofreq = [fA fC fG fT];
  [totx2,expfreq] = repeatsf(monofreq,seq,obsfreq,monoseq);

  obs_exp_freq = [obsfreq expfreq]

  close all;
  i = (ff < 5000);
  plot(f(i),ff(i));
  putbnd(f(i),ff(i));
  putxlab('Frequency of nucleotide G');
  putylab('Total chi-square deviation');

  monoseq = char(monoseq);    % Convert unique identifiers to ascii characters

  return;



