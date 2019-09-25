% ANOVANSTF: Objective function for anovanst: nested unbalanced anova for a 
%           single variable.  
%
%   Syntax: [F,pr,df,ss,ms,varcomp,varprop] = ...
%             anovanstf(X,grps,subgrps,G,SG,n4,getpr,getvc)
%
%         X =       [n x p] matrix of observations.
%         grps =    [n x 1] vector of group labels.
%         subgrps = [n x 1] vector of subgroup labels.
%         G =       [n x ngrps] design matrix for groups.
%         SG =      [n x nsubgrps] design matrix for subgroups.
%         n4 =      within-group sample size.
%         getpr =   boolean variable indicating that significance levels are to 
%                     be estimated [default = 0].
%         getvc =   boolean variable indicating that variance components are 
%                     to be estimated [default = 0].
%         -----------------------------------------------------------------------
%         F =       [2 x p] matrix of F-statistics (among-group, among-subgroup).
%         pr =      [2 x p] matrix of corresponding significance levels, 
%                     asymptotic (if iter=0) or randomized (if iter>0).
%         df =      [4 x p] matrix of degrees of freedom (among-group, 
%                     among-subgroup, within-subgroup, total).
%         ss =      [4 x p] matrix of sums-of-squares (among-group, within-group,
%                     total).
%         ms =      [3 x p] matrix of mean-squares (among-group, 
%                     among-subgroup, within-subgroup).
%         varcomp = [3 x p] matrix of variance-component estimates (among-group,
%                     among-subgroup, within-subgroup).
%         varprop = [3 x p] matrix of variance components as proportions of
%                     total (among-group, among-subgroup, within-subgroup).
%

% Sokal & Rohlf, 1981, pp. 294-299.

% RE Strauss, 1/22/00 - modified from anovanst.
%   2/6/00 - do anovas for series of variables rather than just one.
%   12/3/00 - correct division-by-zeros problem with msS.
%   2/21/02 - check for dfE==0 when calculating msE.

function [F,pr,df,ss,ms,varcomp,varprop] = ...
            anovanstf(X,grps,subgrps,G,SG,n4,getpr,getvc)
  ngrps = size(G,2);
  [n,p] = size(X);

  F = zeros(2,p);
  pr = zeros(2,p);
  df = zeros(4,p);
  ss = zeros(4,p);
  ms = zeros(3,p);
  varcomp = zeros(3,p);
  varprop = zeros(3,p);

  prev_miss = 0;
  G_save = G;
  SG_save = SG;
  grps_save = grps;
  subgrps_save = subgrps;
  n4_save = n4;

  ngrps = size(G,2);
  nsubgrps = size(SG,2);

  for v = 1:p                           % Cycle thru variables
    x = X(:,v);                           % Isolate current variable

    if (prev_miss)                        % If previously had missing data,
      prev_miss = 0;                      %   restore original values
      G = G_save;
      SG = SG_save;
      grps = grps_save;
      subgrps = subgrps_save;
      n4 = n4_save;
      n = length(x);
      ngrps = size(G,2);
      nsubgrps = size(SG,2);
    end;

    i = find(~isfinite(x));               % Delete any missing data
    if (~isempty(i))
      x(i) = [];
      G(i) = [];
      grps(i) = [];
      subgrps(i) = [];

      G = design(grps);                     % Get new group design matrix

      id = 0;                               % Recalc n4
      newsubgrps = zeros(size(subgrps));
      n4 = 0;
      for g = 1:ngrps
        gobs = find(G(:,g));
        sg = uniquef(subgrps(gobs),1);
        n4p = 0;
        for i = 1:length(sg)
          id = id+1;
          cursg = gobs(find(subgrps(gobs)==sg(i)));
          lencur = length(cursg);
          newsubgrps(cursg) = id*ones(lencur,1);
          n4p = n4p + lencur*lencur;
        end;
        n4 = n4 + n4p/length(gobs);
      end;
      subgrps = newsubgrps;

      SG = design(subgrps);                 % Get new subgrp design matrix

      ngrps = size(G,2);
      nsubgrps = size(SG,2);
    end;

    % ANOVA

    totmean = mean(x)*ones(n,1);        % Matching vector of grand means
    y = x - totmean;                    % Deviations from grand mean
    ssT = y'*y;                         % Total sum-of-squares

    gbar = inv(G'*G)*G'*x;              % Group means
    xbar = G*gbar;                      % Matching vector of group means
    e = x - xbar;                       % Deviations from group means
    ssG = e'*e;                         % Within-group sum-of-squares

    gbar = inv(SG'*SG)*SG'*x;           % Subgroup means
    xbar = SG*gbar;                     % Matching vector of subgroup means
    e = x - xbar;                       % Deviations from subgroup means
    ssE = e'*e;                         % Within-subgroup sum-of-squares
    ssS = ssG-ssE;                      % Subgroups within groups sum-of-squares

    ssG = ssT - ssS - ssE;              % Among-subgroup sum-of-squares

    dfT = n-1;                          % Total df
    dfG = ngrps-1;                      % Among-group df
    dfS = nsubgrps-ngrps;               % Among-subgroup df
    dfE = dfT - dfG - dfS;              % Within-subgroup df

    msG = ssG / dfG;
    msS = ssS / dfS;
    if (dfE > 0)
      msE = ssE / dfE;
    else
      msE = NaN;
    end;

    ss(:,v) = [ssG; ssS; ssE; ssT];
    df(:,v) = [dfG; dfS; dfE; dfT];
    ms(:,v) = [msG; msS; msE];

    % Variance components

    ng = diag(G'*G);                      % Group sample sizes
    nsg = diag(SG'*SG);                   % Subgroup sample sizes

    n1 = dfT + 1;
    n2 = nsg'*nsg;
    n3 = ng'*ng;

    np0 = (n4-(n2/n1))/dfG;
    n0 =  (n1-n4)/dfS;
    nb0 = (n1-(n3/n1))/dfG;

    if (getvc)
      varE = msE;
      varS = (msS-msE)/n0;
      varG = (msG-msE-(np0*varS))/nb0;
      varT = varE + varS + varG;

      varcomp(:,v) = abs([varG; varS; varE]);
      varT = sum(varcomp(:,v));
      varprop(:,v) = abs([varG/varT; varS/varT; varE/varT]);
    end;

    % F-statistics

    R = 0;
    C = 0;
    if (np0 ~= n0)
      R = (np0/(np0-n0)) * (msS/msE);
      C = finv(0.975,dfE,dfS) * finv(0.5,dfS,dfE);
    end;

    if (dfS<100 & dfS<(2*dfE) & R>C)    % Use Satterthwaite approximation
      w2 = np0/n0;
      w1 = 1-w2;
      msSp = w1*msE + w2*msS;
      dfSp = (msSp*msSp) / ((w1*w1*msE*msE/dfE) + (w2*w2*msS*msS/dfS));
      Fg = msG / msSp;
      v2 = [floor(dfSp); 0];
    else                                %   else use simple approximate test
      if (msS>0)
        Fg = msG / msS;
      else
        Fg = NaN;
      end;
      v2 = [dfS; 0];
    end;

    Fs = msS / msE;
    F(:,v) = [Fg; Fs];

    % Significance levels

    if (getpr)
      v1 = [dfG; dfS];
      v2(2) = dfE;

      pr(:,v) = 1 - fcdf(F(:,v),v1,v2);
    end;
  end;

  return;
