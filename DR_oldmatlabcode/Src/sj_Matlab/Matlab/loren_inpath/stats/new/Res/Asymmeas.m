% ASYMMEAS: Palmer-Strobeck mixed-model balanced anova model for assessing 
%           bilateral asymmetry, given repeated measures on the same individuals
%           to assess measurement error.
%           (Palmer & Strobeck 1986, ARES 17:405)
%
%   Syntax: [F,pr,df,ss,ms,varcomp,varprop] = asymmeas(x,individual,side)
%
%         x =          [n x p] data matrix for n observations and p variables.
%         individual = [n x 1] classification variable (factor 2, random effect).
%         side =       [n x 1] classification variable (factor 1, fixed effect).
%         -----------------------------------------------------------------------
%         F =       [3 x p] F-statistics.
%         pr =      [3 x p] significance levels of the test.
%         df =      [5 x p] degrees of freedom.
%         ss =      [5 x p] sums-of-squares.
%         ms =      [4 x p] mean-squares.
%         varcomp = [4 x p] variance-component estimates.
%         varprop = [3 x p] variance components as proportions of 
%                           (sides + interaction + error).
%
%         Output components:
%                   (1) sides (directional asymmetry)
%                   (2) individuals (size or shape variation)
%                   (3) interaction (nondirectional asymmetry)
%                   (4) measurement error
%                   (5) total
%

% RE Strauss, 4/30/98

function [F,pr,df,ss,ms,varcomp,varprop] = asymmeas(x,individual,side)

  if (nargout > 2)
    resvect = 1;
  else
    resvect = 0;
  end;
  if (nargout > 5)
    varcomps = 1;
  else
    varcomps = 0;
  end;

  [ntot,p] = size(x);
  if (ntot==1)                        % Input must be column vector
    x = x';
    [ntot,p] = size(x);
  end;

  F = zeros(3,p);
  pr = zeros(3,p);
  df = zeros(5,p);
  ss = zeros(5,p);
  ms = zeros(4,p);
  varcomp = zeros(4,p);
  varprop = zeros(3,p);

  key = individual*(max(side)+1) + side;   % Sort by sides within individuals
  [key,i] = sort(key);
  x = x(i,:);
  side = side(i);
  individual = individual(i);

  [sideval,na] = uniquef(side);        % Number of sides (factor A)
  [indval,nb] =  uniquef(individual);  % Number of individuals (factor B)
  [cellval,n] = uniquef(key);          % Number of sides*individuals (interaction)

  nside = length(sideval);
  nind =  length(indval);
  ncell =  length(cellval);

  if (nside > 2)
    error('  More than two sides represented.');
  end;

  if (any((n-mean(n))~=0))
    disp('  Warning: unbalanced design');
  end;

  n = mean(n);
  na = mean(na);
  nb = mean(nb);

  dfto = ntot-1;                        % Degrees of freedom
  dfa  = nside - 1;
  dfb  = nind - 1;
  dfab = dfa*dfb;
  dfe  = dfto - dfa - dfb - dfab;

  % Statistics by variable (Sokal & Rohlf 1981:333-337; Zar 1996:244-248)

  for pp = 1:p
    xx = x(:,pp);

    sidesum = zeros(nside,1);
    for a = 1:nside
      sidesum(a) = sum(xx(side==sideval(a)));
    end;

    indsum = zeros(nind,1);
    for b = 1:nind
      indsum(b) = sum(xx(individual==indval(b)));
    end;

    cellsum = zeros(ncell,1);
    for ab = 1:ncell
      cellsum(ab) = sum(xx(key==cellval(ab)));
    end; 

    totmean = mean(xx);
    totdev = xx - totmean;

    cellmean = cellsum/n;
    celldev =  cellmean - totmean;

    sidemean = sidesum/na;
    sidedev =  sidemean - totmean;

    indmean = indsum/nb;
    inddev =  indmean - totmean;

    ssto =   (totdev'*totdev);
    sscell = (celldev'*celldev) * n;
    ssa =    (sidedev'*sidedev) * na;
    ssb =    (inddev'*inddev) * nb;
    ssab =   sscell - ssa - ssb;
    sse =    ssto - sscell;

    msa  = ssa / dfa;
    msb  = ssb / dfb;
    msab = ssab / dfab;
    mse  = sse / dfe;

    Fa =  msa / msab;                     
    Fb =  msb / msab;
    Fab = msab / mse;

    % Significance level of observed F-statistics

    f =  [Fa Fb Fab];
    v1 = [dfa dfb dfab];
    v2 = [dfab dfe dfe];
    prob = 1 - fcdf(f,v1,v2);

    F(:,pp) = f';
    pr(:,pp) = prob';

    % Estimates of variance components

    if (varcomps)
      vare =  mse;
      varab = abs(msab-vare)/n;
      varb =  abs(msb-msab)/nb;
      vara =  abs(msa-mse-n*varab)/na;
      varto = vara + varab + vare;

      varcomp(:,pp) = [vara; varb; varab; vare];
      varprop(:,pp) = [vara; varab; vare]/varto;
    end;

    % Output matrices

    if (resvect)
      df(:,pp) = [dfa; dfb; dfab; dfe; dfto];
      ss(:,pp) = [ssa; ssb; ssab; sse; ssto];
      ms(:,pp) = [msa; msb; msab; mse];
    end;

  end;  % next variable

  return;
