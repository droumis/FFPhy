% Reassign: Given two vectors of group labels, reassign the labels of the 
%           second vector to maximize identical matches with the first.  
%           Vectors must be same length.  Numbers of unique labels in the 
%           two vectors can be different.  Solution is not necessarily unique.
%
%     Syntax: [v1,v2] = reassign(v1,v2)
%

% RE Strauss, 12/31/96

function [v1,max_v2] = reassign(v1,v2)
  if (min(size(v1))>1 | min(size(v1))>1)
    error('  Input objects must be vectors.');
  end;
  if (length(v1) ~= length(v2))
    error('  Input vectors must be same length.');
  end;

  v2_init = v2;
  if (size(v2,1)>1)
    colvect = 1;
  else
    colvect = 0;
  end;

  g2 = uniquef(v2);                 % Unique labels in vector 2
  len_g2 = length(g2);

  g2perm = g2;                      % Initialize g2 permutation 
  max_match = sum(v1==v2);          % Initialize matches
  max_v2 = v2;

  for p = 1:(prod(1:len_g2)-1)      % Cycle thru all possible permutations
    g2perm = permnext(g2perm);        % Next permutation

    for i = 1:len_g2                  % Map original vector 2 onto permuted 2
      indx = (v2_init == g2(i));
      if (colvect)
        v2(indx) = g2perm(i)*ones(sum(indx),1);
      else
        v2(indx) = g2perm(i)*ones(1,sum(indx));
      end;
    end;

    match = sum(v1==v2);              % Number of matches
    if (match > max_match)
      max_match = match;
      max_v2 = v2;
    end;
  end;

  return;
