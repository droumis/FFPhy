% ADDMISS:  Add specimens and characters with missing data via a stepwise 
%           procedure so as maintain maximum matrix condition.  Estimates 
%           missing data using the EM algorithm.
%
%     Usage: [Y,obs,vars,propmiss,initsize,initcond,item_added,miss_added,cond_added] = ...
%             addmiss(X,maxmiss,{addchar},{condtype},{matchorig},{varlist},{usecorr})
%
%         X =         [n x p] data matrix having missing data.
%         maxmiss =   maximum proportion of missing data to be estimated.
%         addchar =   optional boolean flag indicating, if true, that variables 
%                       are to be added as well as observations [default = 0].
%         condtype =  optional type of condition factor to be used:
%                       0 = reciprocal condition factor, which measures 
%                           stability under matrix inversion, and is greatest 
%                           when the ratio (eigenvalue 1)/(eigenvalue 2) is  
%                           small [default];
%                       1 = condition factor, which measures stability of the 
%                           eigensolution, and is greatest when the ratio 
%                           (eigenvalue 1)/(eigenvalue 2) is large.
%         matchorig = optional boolean flag indicating, if true, that the 
%                       combinations of observations and variables are to be 
%                       sorted by decreasing deviation from condition factor of 
%                       original matrix [default = 1].
%         varlist =   optional list of variables to use as starting subset (perhaps
%                       a biologically informative list) [default = null].
%         usecorr =   optional boolean flag indicating, if true, that correlation 
%                       matrix is to be used rather than covariance matrix 
%                       [default = 0].
%         ---------------------------------------------------------------------
%         Y =           [m x q] matrix, a subset X for which at most 'maxmiss' data 
%                         have been estimated.
%         obs =         [m x 1] vector of indices of specimens in Y.
%         vars =        [q x 1] vector of indices of characters in Y.
%         propmiss =    proportion of missing data estimated.
%         initsize =    [1 x 2] vector of size of initial complete submatrix.
%         initcond =    scalar indicating the condition of the initial submatrix.
%         item_added =  2-col matrix indicating the index of the item added to 
%                         the initial submatrix; col 1 indicates variables, 
%                         col 2 indicates observations.
%         miss_added =  matching vector indicating the proportion of missing 
%                         data after the item was added to the submatrix.
%         cond_added =  matching vector indicating the matrix condition (rcond)
%                         after the item was added.
%

% RE Strauss, 11/11/00
%   12/18/00 - added optional type of condition factor.
%   6/21/02 -  added option for matching original condition factor.

function [Y,obs,vars,propmiss,initsize,initcond,item_added,miss_added,cond_added] = ...
              addmiss(X,maxmiss,addchar,condtype,matchorig,varlist,usecorr)
  if (nargin < 3) addchar = []; end;
  if (nargin < 4) condtype = []; end;
  if (nargin < 5) matchorig = []; end;
  if (nargin < 6) varlist = []; end;
  if (nargin < 7) usecorr = []; end;

  if (isempty(addchar))
    addchar = 0;
  end;
  if (isempty(condtype))
    condtype = 0;
  end;
  if (isempty(matchorig))
    matchorig = 1;
  end;
  if (isempty(usecorr))
    usecorr = 0;
  end;

  if (maxmiss>1)                        % Convert percentage to proportion
    maxmiss = maxmiss/100;
  end;

  [n,p] = size(X);

  % Get best subset of complete data
%   [nobs,nmiss,combs,condfact] = varcomb(X,condtype); 
  [nobs,nmiss,combvars,combobs,condfact,toomany,C] = ...
                            varcomb(X,condtype,[],[],usecorr,[],matchorig)
                          if (toomany | isempty(combvars))
    error('  ADDMISS: no subset of complete data found.');
  end;

  if (~isempty(varlist))                % Begin with specified list of variables
    c = [combvars; padcols(varlist,combvars)];
    rc = rowtoval(c);
    combpos = find(rc==rc(length(rc)));
    combpos = combpos(1);
    if (isempty(combpos))
      error('  ADDMISS: variable list not found in matrix of combinations.');
    end;
  else
    combpos = 1;
  end;

%   i = find(combvars(combpos,:)>0);
%   vars = combvars(combpos,i)';             % Stash list of vars of subset
%   Y = X(:,vars);                        % Reduce data to vars  
%   r = rowsum(~isfinite(Y));             % Find complete obs
%   obs = find(r==0);
%   Y = Y(obs,:);                         % Reduce data to obs

  obs = combobs(combpos,find(combobs(combpos,:)>0))'
  vars = combvars(combpos,find(combvars(combpos,:)>0))'

  Y = X(obs,vars);                       % Reduce data to complete subset
  initsize = size(Y);
  initcond = condfact(combpos);
  
  if (usecorr)
    initcond = condfactor(C,condtype);
  else
    initcond = condfactor(C,condtype);
  end;

%   if (usecorr)
%     initcond = log(rcond(corrcoef(Y)));
%   else
%     initcond = log(rcond(cov(Y)));
%   end;

  item_added = [];                      % Allocate output matrices
  miss_added = [];
  cond_added = [];

  if (length(obs)==n)
    return;
  end;

  nvars = length(vars);
  nobs = length(obs);

  propmiss = 0;
  totnmiss = 0;
  lowval = -1e6;

  % Add only observations, holding number of variables constant

  if (~addchar)                         
    while (propmiss <= maxmiss)           % Continue till reach max prop missing data
      toadd = 1:n;                          % List of obs to be added
      toadd(obs) = [];
      ntoadd = length(toadd);
      obscond = lowval*ones(ntoadd,1);
      nmiss = zeros(ntoadd,1);
      Xobs = zeros(ntoadd,nvars);

      for i = 1:ntoadd                      % Cycle thru obs to be added
        Xi = X(toadd(i),vars);                % Get current obs
        if (any(isfinite(Xi)))                % If not all missing,
          Yi = [Y; Xi];                         % Append to subset
          Yi = missem(Yi);                      % Estimate missing values
          if (usecorr)
            C = corrcoef(Yi);
          else
            C = cov(Yi);
          end;
          obscond(i) = condfactor(C,condtype);  % Stash cond factor for current obs
          if (matchorig)
            obscond(i) = -abs(obscond(i)-initcond);
          end;
%           obscond(i) = log(rcond(C));           
          Xobs(i,:) = Yi(nobs+1,:);             % And its complete values
          nmiss(i) = sum(~isfinite(Xi));        % Stash amount of missing data
        end;
      end;

      [m,i] = max(obscond);                   % Find obs giving max condition
      if (m > lowval)
        nobs = nobs+1;
        totnmiss = totnmiss + nmiss(i);         % Tally total amount of missing data
        propmiss = totnmiss/(nvars*nobs);
      else
        propmiss = maxmiss+1;
      end;

      if (propmiss <= maxmiss)
        obs = [obs; toadd(i)];                  % Add obs to list
        Y = [Y; Xobs(i,:)];                     % Append to Y

        item_added = [item_added; 0 toadd(i)];
        miss_added = [miss_added; propmiss];
        cond_added = [cond_added; m];
      end;
    end;  % while

    obs = sort(obs);                        % Sort list of obs & vars
    vars = sort(vars);
    Y = X(obs,vars);                        % Restore missing data
    propmiss = sum(sum(~isfinite(Y)))/(length(vars)*length(obs)); % Proportion of missing data
    Y = missem(Y);                          % Predict missing values as a set
  end;

  % Add variables and observations

  if (addchar)
    while (propmiss <= maxmiss)           % Continue till reach max prop missing data
      vartoadd = 1:p;                       % List of vars to be added
      vartoadd(vars) = [];
      ptoadd = length(vartoadd);
      varcond = lowval*ones(ptoadd,1);
      nmiss = zeros(ptoadd,1);
      Xvar = zeros(nobs,ptoadd);

      for ip = 1:ptoadd                     % Cycle thru vars to be added
        Xi = X(obs,vartoadd(ip));            % Get current var
        if (any(isfinite(Xi)))                % If not all missing,
          Yi = [Y Xi];                          % Append to subset
          Yi = missem(Yi);                      % Estimate missing values
          if (usecorr)
            C = corrcoef(Yi);
          else
            C = cov(Yi);
          end;
          varcond(ip) = condfactor(C,condtype);  % Stash cond factor for current obs
          if (matchorig)
            varcond(ip) = -abs(varcond(ip)-initcond);
          end;
%           varcond(ip) = log(rcond(C));          % Stash cond factor for current obs
          Xvar(:,ip) = Yi(:,nvars+1);           % And its complete values
          nmiss(ip) = sum(~isfinite(Xi));       % Stash amount of missing data
        end;
      end;

      [m,i] = max(varcond);                   % Find obs giving max condition
      if (m > lowval)
        nvars = nvars+1;
        totnmiss = totnmiss + nmiss(i);         % Tally total amount of missing data
        propmiss = totnmiss/(nvars*nobs);
      else
        propmiss = maxmiss+1;
      end;

      if (propmiss <= maxmiss)
        vars = [vars; vartoadd(i)];             % Add var to list
        Y = [Y Xvar(:,i)];                      % Append to Y

        item_added = [item_added; vartoadd(i) 0];
        miss_added = [miss_added; propmiss];
        cond_added = [cond_added; m];
      end;

      rca = size(cond_added,1);
      if (rca > 1)
        last_cond = cond_added(rca-1);
      else
        last_cond = initcond;
      end;

      addobs = 0;
      if (propmiss < maxmiss)
        addobs = 1;
      end;
      nadded = 0;

      while (addobs)                        % Add specimens to compensate
        toadd = 1:n;                          % List of obs to be added
        toadd(obs) = [];
        ntoadd = length(toadd);
        obscond = lowval*ones(ntoadd,1);
        nmiss = zeros(ntoadd,1);
        Xobs = zeros(ntoadd,nvars);

        for i = 1:ntoadd                      % Cycle thru obs to be added
          Xi = X(toadd(i),vars);                % Get current obs
          if (any(isfinite(Xi)))                % If not all missing,
            Yi = [Y; Xi];                         % Append to subset
            Yi = missem(Yi);                      % Estimate missing values
            if (usecorr)
              C = corrcoef(Yi);
            else
              C = cov(Yi);
            end;
            obscond(i) = log(rcond(C));           % Stash cond factor for current obs
% Xobs
% Yi
% i
% nobs

%             Xobs(i,:) = Yi(nobs+1,:);             % And its complete values
            Xobs(i,:) = Yi(end,:);                % And its complete values
            nmiss(i) = sum(~isfinite(Xi));        % Stash amount of missing data
          end;
        end;

        [m,i] = max(obscond);                   % Find obs giving max condition
        if (m > lowval)
          nobs = nobs+1;
          totnmiss = totnmiss + nmiss(i);         % Tally total amount of missing data
          propmiss = totnmiss/(nvars*nobs);
        else
          propmiss = maxmiss+1;
          addobs = 0;
        end;

        if (propmiss <= maxmiss)
          obs = [obs; toadd(i)];                  % Add obs to list
          Y = [Y; Xobs(i,:)];                     % Append to Y
          nadded = nadded+1;

          item_added = [item_added; 0 toadd(i)];
          miss_added = [miss_added; propmiss];
          cond_added = [cond_added; m];

          if (nadded >= ceil(n/p) | (m >= last_cond))
            addobs = 0;
          end;
        end;
      end; % while
    end; % while

    vars = sort(vars);                      % Sort lists of vars and observations
    obs = sort(obs);                    
    Y = X(obs,vars);                        % Restore missing data
    propmiss = sum(sum(~isfinite(Y)))/(nvars*nobs); % Proportion of missing data
    Y = missem(Y);                          % Predict missing values as a set
  end;
  
  cond_added = abs(cond_added) + initcond;

  return;
