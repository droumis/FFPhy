% MISSPCA:  Estimates missing data by multiple regression on principal components
%           of complete data.  It is assumed that input data have already been 
%           log-transformed, if appropriate.
%
%     Usage: [Y,propmiss,failflag] = misspca(X,{usecorr})
%
%         X =         [n x p] matrix of original values, with missing data 
%                       indicated by non-finite values (NaN, Inf, -Inf).
%         usecorr =   optional boolean flag indicating, if true, that principal 
%                       components are to be estimated from the correlation 
%                       matrix [default = 0 = covariance matrix].
%         --------------------------------------------------------------------------
%         Y =         [n x p] matrix with missing data replaced by estimated 
%                       values.
%         propmiss =  proportion of missing values in original matrix.
%         failflag =  boolean flag indicating, if true, that the estimation
%                       procedure couldn't be completed due to too much missing data.
%

% RE Strauss, 10/26/00 - modified for Matlab v5 from original function written by 
%                        Joao Alves de Oliveira for Matlab v4.
%

%   12/11/01 - added flag to suppress printing of warnings by pcacov() & pcacorr().
%   2/26/02 -  added failure flag.
%   3/3/02 -   added check for too few complete data.
%   10/1/02 -  modified output argument sequence of linregr().
%   10/8/02 -  modified output argument sequence of linregr().

function [X,propmiss,failflag] = misspca(X,usecorr)
  if (nargin < 2) usecorr = []; end;

  if (isempty(usecorr))
    usecorr = 0;
  end;
  failflag = 0;
  [n,p] = size(X);

  misspos = ~isfinite(X);               % Boolean positions of missing values
  nmiss_obs = rowsum(misspos);          % Numbers of missing values per obs
  u_nmiss_obs = uniquef(nmiss_obs,1);   % Unique numbers of missing vals per obs

  if (sum(nmiss_obs)==0)                % If no missing data, return
    propmiss = 0;
    return;
  else
    propmiss = sum(nmiss_obs)/(n*p);    % Proportion of missing data
  end;

  for iun = 1:length(u_nmiss_obs)       % Cycle thru unique numbers
    if (~failflag)
      icmpl = find(nmiss_obs==0);           % Indices of complete obs
      xcmpl = X(icmpl,:);                   % Isolate submatrix of complete data
      lencmpl = length(icmpl);              % Number of complete obs
      if (lencmpl<3)
        failflag = 1;
      end;
  
      ipred = find(nmiss_obs == u_nmiss_obs(iun));  % Indices for obs to be predicted
      xpred = X(ipred,:);                   % Isolate submatrix of incomplete data
  
      combvar = rowtoval(~isfinite(xpred)); % Identify unique combinations 
      ucombvar = uniquef(combvar);          %   of missing vars
  
      for icv = 1:length(ucombvar)          % Cycle thru combs of missing vars
        if (~failflag)
          i = find(combvar == ucombvar(icv));   % Obs having current combination
          cv = find(isfinite(xpred(i(1),:)));   % List of complete vars
          mv = 1:p;
          mv(cv) = [];                          % List of incomplete vars
  
          xc = [xcmpl(:,cv); xpred(i,cv)];      % Complete submatrix for complete vars
          if (size(xc,1)<3 | rank(cov(xc))<2)
            failflag = 1;
          end;
        
          if (~failflag)
            if (usecorr)                          % PCA of this matrix
              [loadings,percvar,scores] = pcacorr(xc,p,[],[],[],[],1);
            else
              [loadings,percvar,scores] = pcacov(xc,p,[],[],[],[],1);
            end;
      
            scorescmpl = scores(1:lencmpl);             % Scores for complete obs
            scorespred = scores(lencmpl+1:size(xc,1));  % Scores for incomplete obs
            [b,stats,predict] = linregr(scorescmpl,xcmpl(:,mv),[],scorespred);  % Predicted vars
            xpred(i,mv) = predict;                % Plug these predicted values into matrix
          end;
        end;
      end;

      X(ipred,:) = xpred;                     % Plug all predicted values into matrix
      nmiss_obs(ipred) = zeros(length(ipred),1);
    end;
  end;

  return;
