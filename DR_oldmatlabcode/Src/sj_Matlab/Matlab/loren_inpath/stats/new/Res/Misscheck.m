% MISSCHECK: Check input matrices for missing (non-finite) values.
%
%     Usage: [ismiss,matid] = ...
%                 misscheck(p1,{p2},{p3},{p4},{p5},{p6},{p7},{p8},{p9},{p10})
%
%         p1 to p10 - from 1-10 input matrices.
%         --------------------------------------------------------------------
%         ismiss -    boolean value indicating, if true, that at least one 
%                       missing value was found in at least one input matrix.
%         matid -     vector of indices of matrices having one or more missing 
%                       values; null if no missing values found.
%

% RE Strauss, 6/16/00

function [ismiss,matid] = misscheck(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
  ismiss = 0;
  matid = [];

  if (any(~isfinite(p1(:))))
    ismiss = 1;
    matid = [matid; 1];
  end;

  if (nargin >= 2)
    if (any(~isfinite(p2(:))))
      ismiss = 1;
      matid = [matid; 2];
    end;
  end;

  if (nargin >= 3)
    if (any(~isfinite(p3(:))))
      ismiss = 1;
      matid = [matid; 3];
    end;
  end;

  if (nargin >= 4)
    if (any(~isfinite(p4(:))))
      ismiss = 1;
      matid = [matid; 4];
    end;
  end;

  if (nargin >= 5)
    if (any(~isfinite(p5(:))))
      ismiss = 1;
      matid = [matid; 5];
    end;
  end;

  if (nargin >= 6)
    if (any(~isfinite(p6(:))))
      ismiss = 1;
      matid = [matid; 6];
    end;
  end;

  if (nargin >= 7)
    if (any(~isfinite(p7(:))))
      ismiss = 1;
      matid = [matid; 7];
    end;
  end;

  if (nargin >= 8)
    if (any(~isfinite(p8(:))))
      ismiss = 1;
      matid = [matid; 8];
    end;
  end;

  if (nargin >= 9)
    if (any(~isfinite(p9(:))))
      ismiss = 1;
      matid = [matid; 9];
    end;
  end;

  if (nargin >= 10)
    if (any(~isfinite(p10(:))))
      ismiss = 1;
      matid = [matid; 10];
    end;
  end;

  return;
