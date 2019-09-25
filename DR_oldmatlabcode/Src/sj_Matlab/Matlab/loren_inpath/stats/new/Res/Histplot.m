% HISTPLOT: Graphic array of histograms for a set of variables (cols)
%
%     Syntax: histplot(X,{[nrows,ncols]},{fixed_axis})
%
%       X =             [n x p] data matrix.
%       [nrows,ncols] = optional vector of numbers of rows and columns of
%                         histograms to display.  Default = square array (p<50)  
%                         or sets of 10 cols (P>=50).
%       fixed_axis =    optional boolean flag indicating, if true, that scaling 
%                         units are to be preserved for all variables 
%                         [default=0: axes are individually scaled].
%

% RE Strauss, 7/10/98
%   11/25/99 - modified to call histgram().

function histplot(X,nsize,fixed_axis)
  if (nargin < 2) nsize = []; end;
  if (nargin < 3) fixed_axis = []; end;

  [nobs,nvars] = size(X);

  if (isempty(fixed_axis))
    fixed_axis = 0;
  end;

  if (isempty(nsize))                 % Use default array dimensions
    if (nvars < 50)                   % Default is square array of subplots
      ncols = ceil(sqrt(nvars));      %   for fewer than 50 vars,
      nrows = ceil(nvars/ncols);
    else                              % Or sets of 10 columns
      ncols = 10;                     %   for 50 or more vars
      nrows = ceil(nvars/10);
    end;
  else                              % Else check for valid given dimensions
    nrows = nsize(1);
    ncols = nsize(2);
    if ((nrows.*ncols) < nvars)       % Increase nrows if array too small
      nrows = ceil(nvars./ncols);
    end;
  end;

  if (fixed_axis)
    X = X(:);
    g = makegrps(1:nvars,nobs);
    histgram(X,g,tostr(1:nvars),[],[nrows ncols]);
  else
    k = 0;
    for i = 1:nrows
      for j = 1:ncols
        k = k+1;
        if (k<=nvars)
          subplot(nrows,ncols,k);
          histgram(X(:,k),[],[],[],[],[],[],[],[],1);
          puttext(0,1.1,tostr(k),[],12);
        end;
      end;
    end;
  end;

  return;
