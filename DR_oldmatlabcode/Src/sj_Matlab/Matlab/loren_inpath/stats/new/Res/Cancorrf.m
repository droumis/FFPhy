% CANCORRF: Details of canonical correlation analysis.
%
%     Syntax: [r,x_load,y_load,x_scores,y_scores] = cancorrf(X,Y,ncorr,loadflag)
%
%         X =        [n x p] matrix for first suite of variables.
%         Y =        [n x q] matrix for second suite of variables.
%         ncorr =    number of canonical variates [default = min(p,q)].
%         loadtype = optional boolean flag indicating the scaling for the 
%                      loadings: 
%                        0: vector correlations [default];
%                        1: regression coefficients;
%                        2: scoring coefficients.
%         -------------------------------------------------------------------------
%         r =        [c x 1] vector of canonical correlations 
%                      [where c = min(p,q)].
%         x_load =   [p x ncorr] matrix of X variables (rows) with structural 
%                      correlations of canonical variates (cols) 
%                      for first set of variables.
%         y_load =   [q x ncorr] matrix of Y variables (rows) with structural 
%                      correlations of canonical variates (cols)
%                      for second set of variables.
%         x_scores = [n x ncorr] matrix of canonical-variate scores for first set 
%                      of variables.
%         y_scores = [n x ncorr] matrix of canonical-variate scores for second set 
%                      of variables.
%

% RE Strauss, 11/7/01 - isolated from cancorr().
%   11/7/01 - added option of loadtype==2.
%   3/18/02 - standardize scores only if loadtype==0.

function [r,x_load,y_load,x_scr,y_scr] = cancorrf(X,Y,ncorr,loadtype)
  r_only = 1;
  if (nargout > 1)
    r_only = 0;
  end;

  [nobs,p] = size(X);
  [m,q] = size(Y);
  ncvar = q;                              % Number of pairs of canonical variables

  A = corrcoef(X);                        % Correlation matrices
  B = corrcoef(Y);
  C = corr(X,Y);

  Ainv = inv(A);
  R = inv(B)*C'*Ainv*C;
  [b,evals] = eigen(R);                   % Eigenanalysis of correlation matrix

  r = sqrt(evals);                        % Canonical correlations
  if (ncorr < ncvar)
    r = r(1:ncorr);
  end;

  if (~r_only)
    a = zeros(p,p);
    x_scr = zeros(nobs,ncvar);
    y_scr = zeros(nobs,ncvar);

    for i = 1:ncvar                         % Scores on canonical variables
      a(:,i) = Ainv*C*b(:,i);
      for j = 1:nobs
        x_scr(j,i) = a(:,i)'*X(j,:)';
        y_scr(j,i) = b(:,i)'*Y(j,:)';
      end;
    end;

    if (loadtype==0)                        % If loadings are not scoring coefficients
      x_scr = zscore(x_scr);                  % Standardize scores on each can var
      y_scr = zscore(y_scr);
    end;

    switch(loadtype)                        % Loadings
      case 0,                                 % Vector correlations
        x_load = corr(x_scr,X)';            
        y_load = corr(y_scr,Y)';
  
      case 1,                                 % Regression coefficients
        x_load = zeros(p,ncvar);
        y_load = zeros(q,ncvar);
        for ic = 1:ncvar
          S = [ones(nobs,1) x_scr(:,ic)];
          for ip = 1:p
            x = X(:,ip);
            bb = inv(S'*S)*S'*x;
            x_load(ip,ic) = bb(2);
          end;
        end;
        for ic = 1:ncvar
          S = [ones(nobs,1) y_scr(:,ic)];
          for iq = 1:q
            y = Y(:,iq);
            bb = inv(S'*S)*S'*y;
            y_load(iq,ic) = bb(2);
          end;
        end;

      case 2,
        x_load = a;                           % Scoring coefficients
        y_load = b;

      otherwise
        error('  CANCORRF: invalid loadings type');
    end;

    if (ncorr < ncvar)
      ncvar = ncorr;
      x_load = x_load(:,1:ncorr);
      y_load = y_load(:,1:ncorr);
      x_scr = x_scr(:,1:ncorr);
      y_scr = y_scr(:,1:ncorr);
    end;

    for i = 1:ncvar                         % Adjust directions of vectors
      rx = 1;
      ry = 1;
      if (p>1)
        rx = corr(x_load(:,i),ones(p,1),2);   % Check vect corr with vector of ones
      end;
      if (q>1)
        ry = corr(y_load(:,i),ones(q,1),2);
      end;
  
      rxy = [rx,ry];
      [maxabsr,im] = max(abs(rxy));
      s = sign(rxy(im));
  
      if (s<0)                                % If sign of largest abs vector is neg,
        x_load(:,i) = -x_load(:,i);           %   reverse both axes
        y_load(:,i) = -y_load(:,i);
        x_scr(:,i) = -x_scr(:,i);
        y_scr(:,i) = -y_scr(:,i);
      end;
    end;
  end;

  return;
