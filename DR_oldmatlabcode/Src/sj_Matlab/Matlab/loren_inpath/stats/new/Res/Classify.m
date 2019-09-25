% CLASSIFY: If "unknown" observations are provided, singly or in groups, 
%           classifies them into the "known" groups.  If only "known" groups are 
%           provided, reclassifies each observation into one of the known groups.
%
%           In either case, classification is based on minimum Mahalanobis 
%           distances, either original or size-adjusted, and is based on the 
%           assumption that both known and unknown groups have identical  
%           covariance structures.  The training set can optionally be bootstrapped 
%           to estimate frequency distributions of classification.
%
%     Usage: [results,perclass] = classify(X,grps,{Xu},{ugrps},{sizeadj},{iter})
%
%           X =       [n x p] data matrix (obs x vars) for the n "known" training 
%                       observations.
%           grps =    [n x 1] corresponding vector of group identifiers for k groups.
%           Xu =      optional [m x p] data matrix for the m "unknown" observations.
%           ugrps =   optional vector of group identifiers for "unknown" 
%                       observations.  If passed, the unknowns are classified as 
%                       entire groups, based on Mahalanobis distances among group 
%                       centroids; if not passed, the unknowns are classified 
%                       individually based on Mahalanobis distances between 
%                       individuals and group centroids.
%           sizeadj = optional boolean flag indicating, if true, that Mahalanobis 
%                       distances are to be based on residuals from a "size-free" 
%                       discriminant analysis [default = 0].
%           iter =    optional number of bootstrap iterations [default = 0].
%           ------------------------------------------------------------------------
%           results = [m x k] matrix, for the m "unknown" groups of observations 
%                       and k training groups, specifying the proportion of times 
%                       (original solution + bootstrap iterations) in which each 
%                       unknown (rows) is classified into one of the k groups 
%                       (cols).
%           perclass =[k x k] matrix of percent correct classifications, for known 
%                       groups, specifying percentage of time that observations 
%                       from known groups (rows) are reclassified into the same 
%                       or different groups (cols).
%

% RE Strauss, 7/13/98
%   11/29/99 - changed calling sequence.
%   1/25/00 -  changed error messages; changed name of unique().
%   5/18/00 -  added calculations for percent correct classifications.
%   9/23/00 -  updated help documentation.
%   8/27/01 -  base percent correct classifications on bootstrapped results if iter>0.

function [results,perclass] = classify(X,grps,Xu,ugrps,sizeadj,iter)
  if (nargin < 3) Xu = []; end;
  if (nargin < 4) ugrps = []; end;
  if (nargin < 5) sizeadj = []; end;
  if (nargin < 6) iter = []; end;

  if (isempty(iter))                    % Default input auguments
    boot_iter = 0;
  end;
  if (isempty(sizeadj))
    sizeadj = 0;
  end;

  if (isempty(Xu))                      % Determine whether data on
    classify_knowns = 1;                %   "unknowns" has been passed
    classify_unknowns = 0;
  else
    [m,q] = size(Xu);
    classify_knowns = 0;
    classify_unknowns = 1;
    if (~isempty(ugrps))
      ugrps = ugrps(:);
    end;
  end;
    
  [n,p] = size(X);
  grps = grps(:);
  tgrpvals = uniquef(grps);              % Training-groups identifiers
  ngrps = length(tgrpvals);

  if (length(grps) ~= n)
    disp('  CLASSIFY: group vector and data matrix incompatible');
    disp('            for training observations');
    error(' ');
  end;

  perclass = [];

  if (classify_knowns)                  % Jackknife and reclassify 'known' obs
    results = classifk(X,grps,tgrpvals,sizeadj);
    binres = results;

    if (iter)                           % Bootstrap
      for it = 1:iter
        Xb = bootsamp(X,grps);
        results = results + classifk(Xb,grps,tgrpvals,sizeadj);
      end;
      results = results ./ (iter+1);
      [m,im] = max(results');
      binres = zeros(size(binres));
      for ir = 1:size(binres,1)
        binres(ir,im(ir)) = 1;
      end;
    end;
    
    perclass = zeros(ngrps,ngrps);
    for i = 1:ngrps
      ig = find(grps == tgrpvals(i));
      lg = length(ig);
      perclass(i,:) = 100*sum(binres(ig,:))/lg;
    end;

  end;

  if (classify_unknowns)
    if (p ~= q)
      disp('  CLASSIFY: numbers of variables for training and unknown observations');
      disp('            must be equal');
      error(' ');
    end;

    len_ugrps = length(ugrps);
    if (len_ugrps > 0)
      if (len_ugrps ~= m)
        disp('  CLASSIFY: group vector and data matrix incompatible');
        disp('            for unknown observations');
        error(' ');
      end;
    else
      ugrps = [1:m]';
    end;

    ntgrps = length(tgrpvals);              % Number of training groups
    maxtgrp = max(tgrpvals);                % Max training-group label

    ugrps = ugrps + maxtgrp;                % Ensure that unknown grp-ids are unique
    ugrpvals = uniquef(ugrps);              % Groups of unknowns

    results = classifu(X,grps,Xu,ugrps,ntgrps,ugrpvals,sizeadj);

    if (iter)                               % Bootstrap the training groups
%level = 0.4;                                  % Level of bootstrap subst: 0-1
%f = max([round((1-level)*n),1]);
%f = min([f,n]);

      for it = 1:iter
%it
        bX = bootsamp(X,grps);
%bX(1:f,:) = X(1:f,:);
        results = results + classifu(bX,grps,Xu,ugrps,ntgrps,ugrpvals,sizeadj);
      end;
      results = results./(iter+1);
    end;
  end;

  return;

