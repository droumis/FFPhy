% SIZEFREP: Parse the solution row vector from SIZEFREF().
%
%     Syntax: [loadings,percvar,scores,fscores,D2,R,wload,wperc,wsize,S]
%                     = sizefrep(solution,outmat,n,p,g,ndf)
%
%        solution =    row vector of concatenated results.  Individual output 
%                        matrices are concatenated by column, then 
%                        transposed.
%        outmat =      boolean vector indicating results to be returned:
%                         1) loadings: (vector correlations)        [p x ndf] 
%                         2) percvar:  percents total variance      [ndf x 1] 
%                         3) scores:   DF scores                    [n x ndf]
%                         4) fscores:  floated DF scores            [m x ndf]
%                         5) D2:       Mahal distances              [g x g]
%                         6) R:        size-invariant residuals     [n x p]
%                         7) wload:    size-vector loadings         [p x 1]
%                         8) wperc:    percent size-vector variance [1 x 1]
%                         9) wsize:    within-group size scores     [n x 1]
%                        10) S:        among-group size scores      [n x 1]
%        outsize =     sizes of matrices corresponding to outmat elements.
%        n,p,g,ndf =   number of observations, variables, groups, and 
%                        discriminant functions.
%           

% RE Strauss, 2/9/97
%   5/2/00 -  return floated scores.

function [loadings,percvar,scores,fscores,D2,R,wload,wperc,wsize,S] ...
             = sizefrep(solution,outmat,outsize,n,p,g,ndf)
  loadings = [];
  percvar = [];
  scores = [];
  fscores = [];
  D2 = [];
  R = [];
  wload = [];
  wperc = [];
  wsize = [];
  S = [];

  j = 0;

  if (outmat(1))        % loadings: [p x ndf]
    i = j+1;
    j = (i-1) + outsize(1);
    loadings = reshape(solution(i:j),p,ndf);
  end;
  
  if (outmat(2))        % percvar:  [ndf x 1]
    i = j+1;
    j = (i-1) + outsize(2);
    percvar = solution(i:j)';
  end;
  
  if (outmat(3))        % scores:   [n x ndf]
    i = j+1;
    j = (i-1) + outsize(3);
    scores = reshape(solution(i:j),n,ndf);
  end;
  
  if (outmat(4))        % fscores:  [n x ndf]
    i = j+1;
    j = (i-1) + outsize(4);
    fscores = reshape(solution(i:j),outsize(4)/ndf,ndf);
  end;
  
  if (outmat(5))        % D2:       [g x g]
    i = j+1;
    j = (i-1) + outsize(5);
    D2 = reshape(solution(i:j),g,g);
  end;
  
  if (outmat(6))        % R:        [n x p]
    i = j+1;
    j = (i-1) + outsize(6);
    R = reshape(solution(i:j),n,p);
  end;
  
  if (outmat(7))        % wload:    [p x 1]
    i = j+1;
    j = (i-1) + outsize(7);
    wload = solution(i:j)';
  end;
  
  if (outmat(8))        % wperc:    [1 x 1]
    i = j+1;
    j = (i-1) + outsize(8);
    wperc = solution(i:j);
  end;
  
  if (outmat(9))        % wsize:    [n x 1]
    i = j+1;
    j = (i-1) + outsize(9);
    wsize = solution(i:j)';
  end;
  
  if (outmat(10))       % S:        [n x 1]
    i = j+1;
    j = (i-1) + outsize(10);
    S = reshape(solution(i:j),n,1);
  end;
  
  return;

