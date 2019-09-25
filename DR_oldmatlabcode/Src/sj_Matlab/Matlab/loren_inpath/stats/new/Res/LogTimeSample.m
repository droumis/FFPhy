% LogTimeSample: Given the first and last days of sampling and the number of desired
%                samples, returns the sequence of sampling days on a logarithmic scale.
%
%     Syntax: days = logtimesample(first,last,N)
%
%         first = scalar indicating first sampling day.
%         last =  scalar indicating last sampling day.
%         N =     number of sampling days.
%         ---------------------------------------------
%         days =  row vector of sampling days.
%

% RE Strauss, 6/22/03

function days = logtimesample(first,last,N)
  days = exp(linspace(log(first),log(last),N));

  return;
  