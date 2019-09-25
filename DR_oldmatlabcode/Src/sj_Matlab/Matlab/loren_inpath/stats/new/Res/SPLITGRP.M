% SPLITGRP: Given a vector 'grp' containing grouping codes for 2-9 groups,
%     and a matrix X containing corresponding data, this procedure splits
%     the data matrix X into 2-9 separate matrices, one per group.  Groups are 
%     identified in the sequence in which they are encountered in 'grp'; corresponding 
%     group labels are given in the column vector 'grp_id'.
%
%     Usage: [grp_id,X1,X2,...,X9] = splitgrp(grps,X)
%
%     Matrices X3-X9 are empty if only 2 groups are present, and so on.
%

% RE Strauss.  Rewritten 9/98.
%   11/29/99 - reversed calling sequence.
%   1/25/00 -  revised error message; changed name of unique().

function [index,X1,X2,X3,X4,X5,X6,X7,X8,X9] = splitgrp(X,grps)
  [rx,cx] = size(X);
  [rg,cg] = size(grps);
  if (rx ~= rg),
    error('   Group and data matrices must have same number of rows');
  end;

  X1=[]; X2=[]; X3=[]; X4=[]; X5=[];   % Allocate output matrices
  X6=[]; X7=[]; X8=[]; X9=[];

  index = uniquef(grps);               % Get set of group labels
  ngrps = length(index);

  if (ngrps > 9)
    error('  SPLITGRP: Max of 9 groups allowed.');
  end;

  i = find(grps==index(1));
  X1 = X(i,:);

  if (ngrps>1 & nargout>2)
    i = find(grps==index(2));
    X2 = X(i,:);
  end;

  if (ngrps>2 & nargout>3)
    i = find(grps==index(3));
    X3 = X(i,:);
  end;

  if (ngrps>3 & nargout>4)
    i = find(grps==index(4));
    X4 = X(i,:);
  end;

  if (ngrps>4 & nargout>5)
    i = find(grps==index(5));
    X5 = X(i,:);
  end;

  if (ngrps>5 & nargout>6)
    i = find(grps==index(6));
    X6 = X(i,:);
  end;

  if (ngrps>6 & nargout>7)
    i = find(grps==index(7));
    X7 = X(i,:);
  end;

  if (ngrps>7 & nargout>8)
    i = find(grps==index(8));
    X8 = X(i,:);
  end;

  if (ngrps>8 & nargout>9)
    i = find(grps==index(9));
    X9 = X(i,:);
  end;


  return;
