% PROBDETECT: Estimates the probability of detecting one or more occurrences of an 
%             event, given the sample size and probability of the event occurring 
%             (i.e., the binomial parameters).
%
%     detect = probdetect(N,P)
%             
%         N =       vector (length n) of sample sizes.
%         P =       vector (length p) of probabilities of occurrence.
%         ------------------------------------------------------------------------
%         detect =  [n x p] matrix of probabilities of detection, for all possible
%                     combinations of N and P.
%

% RE Strauss, 12/5/01

function detect = probdetect(N,P)
  N = N(:);
  P = P(:);
  
  n = length(N);
  p = length(P);
  detect = zeros(n,p);
  
  for i = 1:n
    for j = 1:p
      detect(i,j) = 1-binopdf(0,N(i),P(j));
    end;
  end;

  return;
  