% NMATCH: Finds the probability of finding a match of a word within a string of 
%         random letters.
%
%     Usage: [pr,expstrlen] = nmatch(nlet,wordlen,stringlen,{bothdir})
%
%         nlet =      number of letters in alphabet.
%         wordlen =   number of letters in word.
%         stringlen = number of letters in string.
%         bothdir =   optional boolean flag indicating, if true, that the word 
%                       can be found in either the forward or reverse direction 
%                       (e.g., for DNA sequences) [default = 0].
%         ---------------------------------------------------------------------
%         pr =         probability of a match.
%         expstrlen =  expected length of sequence containing 1 match.
%

% RE Strauss, 2/21/00

function [pr,expstrlen] = nmatch(nlet,wordlen,stringlen,bothdir)
  if (nargin < 4) bothdir = []; end;

  if (isempty(bothdir))
    bothdir = 0;
  end;

  lognum = (stringlen-wordlen)*log(nlet)...
           + sum(log(1:stringlen)) - sum(log(1:wordlen))...
           - sum(log(1:(stringlen-wordlen)));
  logden = stringlen*log(nlet);

  pr = exp(lognum-logden);
  if (bothdir)
    pr = pr/2;
  end;
  expstrlen = 1/pr;

  return;
