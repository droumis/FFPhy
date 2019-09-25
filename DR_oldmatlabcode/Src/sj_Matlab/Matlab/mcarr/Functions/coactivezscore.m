function z = coactivezscore(n1, n2, varargin)
%function zscore = coactivez(n1, n2, varargin)
%
% Computes the probability that the events represented by n1 and n2 were
% independent.
%
% n1 and n2 must both be 1xN where each element represents the number of
% events (e.g. spikes) occurred during each of the N windows (e.g. ripples).
% n1 and n2 can also be logicals indicating whether an event occurred.
%
% the zscore represents deviation from the mean number of coactive events
% expected and is calculated as follows:
%
% z = ((num12 - num1 * num2)/Nwindows) / sqrt(num1 * num2 * (NWindows - n1) *
%      (Nwindows - n2) / (NWindows^2 * (NWindows - 1)))
%
%
% Annabelle's version of the function -- has 4 input arguments
%
% function z = coactivezscore(N, nAB, nA, nB)
% % z = coactivezscore(N, nAB, nA, nB)
% % N is number of events
% % nAB is number of events in which A&B were active
% % nA is # events in which A was active
% % nB is # events in which B was active
% %
% % undefined if one cell never fires in ripples or always fires in all
% % ripples (ie N=nA, N=nB, nA=0 or nB=0), or only 1 ripple  (N=1)
%


if nargin == 2
    %Loren's version has 2 input arguments
    % convert n1 and n2 to logicals
    n1 = logical(n1);
    n2 = logical(n2);
    N = length(n1);

    if (length(n2) ~= N)
        error('n1 and n2 must be the same length')
    end

    num1 = sum(n1);
    num2 = sum(n2);
    num12 = sum(n1 & n2);

    V = num1 * num2 * (N - num1) * (N - num2) / (N^2 * (N-1));
    z = (num12 - num1 * num2 / N) / sqrt(V);

elseif nargin ==4  
    %Annabelle's version (called by diff scripts/functions with 4
    %arguments)
    %transfer names of variables
    N =n1;
    nAB = n2;
    nA = varargin{1};
    nB = varargin{2};

    %calculate z-score
    E = nA*nB/N ;
    sig = ( nA*nB*(N-nA)*(N-nB) )/ (N*N*(N-1));

    z = (nAB - E) / sqrt(sig);
end




% Annabelle's version of the function -- hjas 4 input arguments
%
% function z = coactivezscore(N, nAB, nA, nB)
% % z = coactivezscore(N, nAB, nA, nB)
% % N is number of events
% % nAB is number of events in which A&B were active
% % nA is # events in which A was active
% % nB is # events in which B was active
% %
% % undefined if one cell never fires in ripples or always fires in all
% % ripples (ie N=nA, N=nB, nA=0 or nB=0), or only 1 ripple  (N=1)
%
% E = nA*nB/N ;
% sig = ( nA*nB*(N-nA)*(N-nB) )/ (N*N*(N-1));
%
% z = (nAB - E) / sqrt(sig);
%
% end
