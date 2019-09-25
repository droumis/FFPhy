% MAKEGRPS:  Compose a group-membership column vector, given the group 
%            labels and frequencies.
%
%     Syntax:  grps = makegrps(labels,freqs)
%
%            labels = vector of group labels.
%            freqs  = corresponding vector of sample sizes for each group,
%                       or a scalar if group sample sizes are equal.
%            --------------------------------------------------------------
%            grps =   grouping vector.
%

% RE Strauss, 3/27/97
%   9/22/98 - correct error checks

function grps = makegrps(labels,freqs)
  [rl,cl] = size(labels);
  [rf,cf] = size(freqs);

  if (min([rl,cl]) > 1 | min([rf,cf]) > 1)
    error('MAKEGRPS: labels and frequencies must be vectors.');
  end;

  if (max([rf,cf]) == 1)        % Expand scalar for equal sample sizes
    freqs = freqs * ones(rl,cl);
  end;

  if (length(labels) ~= length(freqs))
    error('MAKEGRPS: label and frequency vectors must be same length');
  end;

  grps = zeros(sum(freqs),1);
  last = 0;
  for i = 1:length(labels)
    first = last+1;
    last = first+freqs(i)-1;
    grps(first:last) = labels(i)*ones(freqs(i),1);
  end;

  return;
