% PATRIST: Given an ancestor function, calculates the step matrix among 
%           terminal taxa.  Given an ancestor function and vector of 
%           corresponding branch lengths, also calculates the patristic-
%           distance matrix among terminal taxa.  The step matrix gives
%           the number of resolved nodes between pairs of taxa after
%           considering polytomies.
%
%     Usage: [step,patr] = patrist(anc,{brlen})
%
%           anc =   an ancestor-function vector, or a 2- or 3-column matrix specifying
%                     descendant labels (col 1), ancestor labels (col 2), and optional
%                     branch lengths (col 3).
%           brlen = optional corresponding vector of corresponding branch (edge) 
%                     lengths.
%           --------------------------------------------------------------------------
%           step = step matrix, measuring the number of nodes between pairs 
%                     of terminal taxa.
%           patr = patristic distances among pairs of terminal taxa.
%

% RE Strauss, 7/14/95
%    9/20/99 - update handling of null input arguments.
%   10/18/02 - add error message for ancestor function lacking 0.
%   10/28/02 - call ancestorfunction() to resolve polytomies;
%              adjust step matrix for zero branch lengths.

function [step,patr] = patrist(anc,brlen)
  if (nargin < 2) brlen = []; end;

  if (isempty(brlen))
    anc = ancestorfunction(anc);
  else
    [anc,brlen] = ancestorfunction(anc,brlen);
  end;
  
  n = (length(anc)+1)/2;              % Number of taxa
  step = zeros(n,n);                  % Allocate step matrix

  if (~isempty(brlen))                % If branch lengths have been passed,
    calc_patr = 1;                    %   set flag &
    patr = step;                      %   allocate patristic-distance matrix
  else
    calc_patr = 0;
    patr = [];
  end;
  
  for t1 = 1:(n-1)                    % For all pairs of taxa,
    anc_list = [];                    %   Make list of ancestors of taxon 1
    brlen_list = [];
    a = t1;
    if (calc_patr)
      dist_list = [];                 %   Corresponding cum branch lengths
      cum_dist = 0;
    end;
    while (anc(a)~=0)
      anc_list = [anc_list anc(a)];
      brlen_list = [brlen_list brlen(a)];      
      if (calc_patr)
        cum_dist = cum_dist + brlen(a);
        dist_list = [dist_list cum_dist];
      end;
      a = anc(a);
    end;

    for t2 = (t1+1):n
      t2_steps = 0;
      cum_dist = 0;
      a = t2;
      k = [];
      while (isempty(k))                  % Find first common ancestor with t1
        if (brlen(a)~=0)
          t2_steps = t2_steps+1;
        end;
        cum_dist = cum_dist + brlen(a);
        k = find(anc_list==anc(a));
        if (isempty(k))
          a = anc(a);
        end;
      end;  % End while
      
      
      t1_steps = sum(brlen_list(1:k)>0);
      steps = t2_steps + t1_steps - 1;
      step(t1,t2) = steps;
      step(t2,t1) = steps;
      if (calc_patr)
        patr(t1,t2) = dist_list(k) + cum_dist;
        patr(t2,t1) = dist_list(k) + cum_dist;
      end;
    end;  % End for t2
  end;  % End for t1

  return;
