function z = probco(n1, n2)
%function zscore = coactivez(n1, n2)
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

