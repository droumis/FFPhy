function out = jointsurprise(N, nAB, nA, nB)
% out = jointsurprise(N, nAB, nA, nB)
% see Pazienti & Grun, J Computtional Neuroscience 2006
% joint surprise gies us a p value reflecting the probability that the number 
% of observed coincident events is higher than expected based on baseline activity.
% The probability is computer with a poisson function
% p = 1 - p(number coincident event <nAB)
%   = 1 - sum(n:0-->nAB) ( (npred^n) * (e^-npred) / n! )
%
% INPUTS
%   N = number of events or bins of time
%   nAB = number of joint events, in which A&B were both active
%   nA = number of events in which A was active
%   nB = number of events in which B was active
%
% OUT = joint surprise
%       if nA or nB = 0, joint surprise = 1
%
% ASinger 4/15/09

pA = nA/N;
pB = nB/N;
npred = pA * pB * N;

%comput sum of p( # coincident events < nAB)
ps= poisscdf([0:nAB],npred); 

if nA==0 || nB==0 % special case: if one cell never fired
    ps = 0;
end

out = 1-ps(end);
end