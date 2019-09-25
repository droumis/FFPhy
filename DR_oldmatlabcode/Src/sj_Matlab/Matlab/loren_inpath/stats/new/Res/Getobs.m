% GETOBS: Given a grouping vector for N observations and a list of groups, extracts the 
%         Ng observations for the specified groups.
%
%     Usage: [subsetobs,nf] = getobs(obs,groups)
%
%         obs =       [N x 1] group-membership vector for N observations.
%         groups =    vector of labels of groups to be extracted.
%         -----------------------------------------------------------------------
%         subsetobs = [Ng x 1] vector of the indices (addresses) of the subset of 
%                       observations in the groups specified by 'groups'.
%         nf =        boolean flag indicating that one or more groups were not 
%                       found.
%

% RE Strauss, 6/11/98
%   4/19/00 - added 'not found' flag.

function [subsetobs,nf] = getobs(obs,groups)
  [r,c] = size(groups);
  if (r==1) 
    groups = groups';
  end;

  nf = 0;
  subsetobs = [];
  for g = groups'
    i = find(obs == g);
    if (~isempty(i))
      subsetobs = [subsetobs; i];
    else
      nf = 1;
    end;
  end;

  return;
