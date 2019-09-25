% RANDMISS: Randomly fills a matrix with missing values, ensuring that no rows nor
%           columns have all missing values.
%
%     [Y,toomany,rn,cn] = randmiss(X,nmiss,{rowmin},{colmin})
%
%           X =       arbitrary input matrix.
%           nmiss =   proportion (if <1) or number (if >=1) of missing values 
%                       to be dispersed.
%           rowmin =  minimum number of non-missing values per row [default = 1].
%           colmin =  minimum number of non-missing values per column [default = 
%                       rowmin].
%           ---------------------------------------------------------------------
%           Y =       resultant matrix.
%           toomany = boolean flag for too many missing data, given rowmin and
%                       colmin.
%           rn =      vector of rows of missing values.
%           cn =      matching vector of cols of missing values.
%

% RE Strauss, 12/18/98
%  3/17/99 -  added systematic attempt to correct for saturated rows or columns, 
%               with 'toomany' flag.
%  3/18/99 -  modified function to call randbinm() to find positions for missing 
%               data.
%  4/17/00 -  return indices to missing values.
%  11/13/00 - allow for number of missing values to be a proportion.
%  3/3/02 -   return a 'toomany' flag for too many missing data, rather than print
%               an error message.

function [X,toomany,rn,cn] = randmiss(X,nmiss,rowmin,colmin)
  if (nargin < 3) rowmin = []; end;
  if (nargin < 4) colmin = []; end;

  if (isempty(rowmin))                  % Default input parameters
    rowmin = 1;
  end;
  if (isempty(colmin))
    colmin = rowmin;
  end;

  [nrows,ncols] = size(X);
  ncells = nrows*ncols;

  if (nmiss < 1 )                       % Convert proportion to number
    nmiss = round(nmiss*ncells);
  end;

  ngood = ncells - nmiss;

  toomany = 0;                          % Return flag for too many missing data
  if (ngood<(nrows*rowmin) | ngood<(ncols*colmin))
    toomany = 1;
    X = [];
    rn = [];
    cn = [];
    return;
  end;

  rowtot = ones(1,nrows)*rowmin;        % Row totals for good values
  ng = ngood - nrows*rowmin;
  if (ng>0)
    rg = makegrps(1:nrows,(ncols-rowmin)*ones(1,nrows));
    rg = rg(randperm(length(rg)));
    for i = 1:ng
      rowtot(rg(i)) = rowtot(rg(i))+1;
    end;
  end;

  coltot = ones(1,ncols)*colmin;        % Column totals for good values
  ng = ngood - ncols*colmin;
  if (ng>0)
    cg = makegrps(1:ncols,(nrows-colmin)*ones(1,ncols));
    cg = cg(randperm(length(cg)));
    for i=1:ng
      coltot(cg(i)) = coltot(cg(i))+1;
    end;
  end;

  randpos = randbinm(rowtot,coltot);    % Generate random binary matrix

  [rn,cn] = find(randpos==0);           % Substitute NaN's for zeros in binary matrix
  for i = 1:nmiss
    X(rn(i),cn(i)) = NaN;
  end;

  return;
