% GroupMeans: given a data vector and group-identification vector, returns the
%             within-group means and standard deviations.
%
%     Usage: [grp_ids,m,s] = groupmeans(x,{g})
%
%         x = [n x p] matix of data values for k groups.
%         g = optional corresponding [n x 1] vector of group identifiers
%               [default = vector of ones].
%         ------------------------------------------------------
%         grp_ids = [k x 1] vector of unique group identifiers.
%         m = [k x p] vector of group means.
%         s = [k x p] vector of group standard deviations.
%

function [grp_ids,m,s] = groupmeans(x,g)
  if (nargin <2)
    g = [];
  end;
  if (nargin < 3)
    opt = [];
  end;

  [n,p] = size(x);

  if (isempty(g))
    g = ones(n,1);
  end;

  if (length(g) ~= n)
    error('  GroupMeans: x and g must have same number of observations.');
  end;
  
  grp_ids = uniquef(g);
  k = length(grp_ids);
  
  m = zeros(k,p);
  s = zeros(k,p);
  
  for current_group = 1:k                     % Cycle through groups
    obs = find(g == grp_ids(current_group));
    
    m(current_group,:) = mean(x(obs,:));
    s(current_group,:) = std(x(obs,:));
  end;
  
  return;
  