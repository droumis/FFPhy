% MakeRepeatSeqs: Produces a column vector of k repeated sequences of the integers
%                 1 through n: [1,2,...,n, 1,2,...,n, ...]'
%
%     Usage: seqs = makerepeatseqs(n,k)
%
%         n = max value within sequence.
%         k = number of sets of sequences.
%         -----------------------------------------
%         seqs = resulting [k*n x 1] column vector.
%

% RE Strauss, 3/14/03

function seqs = makerepeatseqs(n,k)
  s = [1:n]' * ones(1,k);
  seqs = s(:);
  
  return;
  