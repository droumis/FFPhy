% SIZEFREF: Objective function for size-invariant discriminant analysis,
%           called by SIZEFREE() and BOOTSTRP().
%
%     Syntax: solution = sizefref(X,grps,Xf,initsoln,nu,ngrps,...
%                                 ndf,within,outmat,outsize,loadtype,kindregr)
%
%        X =           [n x p] data matrix (obs x vars).
%        grps =        row or column vector of group identifiers.
%        Xf =          optional [m x p] matrix of observations to be "floated" 
%                        onto the discriminant functions.
%        initsoln =    initial solution for bootstrapping; used to adjust 
%                      signs of loadings to be consistent with original.
%        nu =          passed by bootstrp() but not used.
%        ngrps =       number of groups.
%        ndf =         optional number of leading discriminant functions for
%                        which scores are desired (default = groups-1).
%        within =      kind of size vector: 1 = within-group, 0 = among-group.
%        outmat =      boolean vector indicating results to be returned:
%                        1) loadings: (vector correlations)        [p x ndf] 
%                        2) percvar:  percents total variance      [ndf x 1] 
%                        3) scores:   DF scores                    [n x ndf]
%                        4) D2:       Mahal distances of residuals [g x g]
%                        5) R:        size-invariant residuals     [n x p]
%                        6) wload:    size-vector loadings         [p x 1]
%                        7) wperc:    percent size-vector variance [1 x 1]
%                        8) wsize:    within-group size scores     [n x 1]
%                        9) S:        among-group size scores      [n x 1]
%        outsize =     sizes of matrices corresponding to outmat elements.
%        loadtype =    optional boolean flag indicating the scaling for the 
%                        loadings: 
%                          0: vector correlations [default];
%                          1: regression)coefficients;
%                          2: squared loadings sum to unity.
%        kindregr -    regression model used: 0 = major axis [default], 1 = predictive.
%        ------------------------------------------------------------------------------
%        solution =    row vector of concatenated results.  Individual output 
%                        matrices are concatenated by column, then 
%                        transposed.
%           

% RE Strauss, 12/9/97
%   5/2/00 -  return floated scores.

function solution = sizefref(X,grps,Xf,initsoln,nu,ngrps,...
                             ndf,within,outmat,outsize,loadtype,kindregr)

  [nobs,nvars] = size(X);       % Numbers of observations & variables

  % Size-invariant residuals from pooled within-group size vector
  if (within)
     Z = grpcentr(X,grps);      % Group-center all variables
  else
     Z = X;                     %   or not
  end;

  [wload,wperc,wsize] = pcacov(Z,1); % Within-group size scores
  [l,p,S,Sf] = pcacov(X,1,0,Xf);     % Among-group size scores


  if (kindregr == 0)            % Regr slopes on within-grp size vector
    B = zeros(2,nvars);
    for i = 1:nvars
      b = majaxis(wsize,Z(:,i));  % Major axis 
      B(1,i) = b(1,2);
      B(2,i) = b(1,1);
    end;
  else
    B = linregr(wsize,Z);         % Predictive regr
  end;

  pred = S * B(2,:);            % Pred values based on among-grp slopes
  X = X-pred;                   % Centered size-invariant residuals
  mX = mean(X);
  R = X - ones(nobs,1)*mX;

  if (~isempty(Xf))             % Same for observations to be floated
    Xf = Xf - Sf*B(2,:);
    Rf = Xf - ones(size(Xf,1),1)*mX;
  else
    Rf = [];
  end;

  % Discriminant analysis on residuals
  loadings = discrim(R,grps,ndf);   % All vars
  [sortl,index] = sort(abs(loadings));
  omit_var = index(1);
  R_omit = R(:,[omit_var]);     % Delete var with least discrim value
  R(:,[omit_var]) = [];         %   to avoid singularity

  if (~isempty(Rf))             % Delete var from floated obs also
    Rf(:,[omit_var]) = [];
  end;

  [loadings,percvar,scores,fscores] = discrim(R,grps,ndf,loadtype,Rf);   % Rerun discrim
  load_omit = corr(R_omit,scores);   % Loading for omitted variable

  if (outmat(5))
    D2 = mahal(R,grps);         % Get Mahalanobis distance
  end;

  scores = real(scores);        % If complex, save real portions
  loadings = real(loadings);
  load_omit = real(load_omit);

  if (omit_var == 1)            % Restore missing loading
    loadings = [load_omit; loadings];
    R = [R_omit, R];
  elseif (omit_var == nvars)
    loadings = [loadings; load_omit];
    R = [R, R_omit];
  else
    loadings = [loadings(1:(omit_var-1),:);
                load_omit;
                loadings(omit_var:(nvars-1),:)];
    R = [R(:,1:(omit_var-1)), ...
         R_omit, ...
         R(:,omit_var:(nvars-1))];
  end;

  % Concatenate results

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

  i = find(outmat);
  solution = zeros(1,sum(outsize(i)));
  cumsize = cumsum(outsize);

  if (outmat(1))        % loadings: loadings as vector correlations
    if (~isempty(initsoln))     % Check signs of loadings
      is = reshape(initsoln(1:cumsize(1)),nvars,ndf);
      for i=1:ndf
        if (corr(loadings(:,i),is(:,i))<0)
          loadings(:,i) = -loadings(:,i);
        end;
      end;
    end;
    solution(1:cumsize(1)) = loadings(:)';
  end;
  
  if (outmat(2))        % percvar:  percentages of total variance
    solution((cumsize(1)+1):cumsize(2)) = percvar';
  end;
  
  if (outmat(3))        % scores:   DF scores
    solution(cumsize(2)+1:cumsize(3)) = scores(:)';
  end;
  
  if (outmat(4))        % fscores:  floated DF scores
    solution(cumsize(3)+1:cumsize(4)) = fscores(:)';
  end;
  
  if (outmat(5))        % D2:       Mahalanobis distances
    solution(cumsize(4)+1:cumsize(5)) = D2(:)';
  end;
  
  if (outmat(6))        % R:        size-invariant residuals
    solution(cumsize(5)+1:cumsize(6)) = R(:)';
  end;
  
  if (outmat(7))        % wload:    size-vector loadings
    if (initsoln~=[])     % Check signs of loadings
      is = initsoln(cumsize(6)+1:cumsize(7))';
      if (corr(wload,is)<0)
        wload = -wload;
      end;
    end;
    solution(cumsize(6)+1:cumsize(7)) = wload';
  end;
  
  if (outmat(8))        % wperc:    size-vector percent total variance
    solution(cumsize(7)+1:cumsize(8)) = wperc';
  end;
  
  if (outmat(9))        % wsize:    within-group size scores
    solution(cumsize(8)+1:cumsize(9)) = wsize';
  end;
  
  if (outmat(10))       % S:        among-group size scores
    solution(cumsize(9)+1:cumsize(10)) = S';
  end;
  
  return;

