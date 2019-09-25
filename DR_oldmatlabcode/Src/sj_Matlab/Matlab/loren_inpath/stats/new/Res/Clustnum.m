% CLUSTNUM: Evaluate putative cluster assignments of points based on ADR, the 
%           log-ratio of mean among-group distance to mean within-group distance.  
%           Only distances along the Gabriel graph are used to calculate distances.
%
%     Usage: [adr,K,tprob,Uprob,df,t,U] = clustnum(crds,grps)
%
%           crds =  [n x 2] matrix of point coordinates.
%           grps =  [n x k] matrix of group-identifier vectors.
%           -----------------------------------------------------------------
%           adr =   [1 x k] vector of ADR statistic values.
%           K =     corresponding vector of number of groups.
%           tprob = corresponding vector of t-probabilities (bias-corrected).
%           Uprob = corresponding vector of U-probabilities, estimated 
%                     as t-probabilities on the rank distances.
%           df =    corresponding vector of degrees of freedom.
%           t =     corresponding vector of t statistics.
%           U =     corresponding vector of U statistics.
%

% RE Strauss, 1/11/99
%   6/11/99 -  modified to output probabilities.
%   10/23/02 - added bias-correction for t probability;
%              return K;
%              make probability calculations optional.

function [adr,K,tprob,Uprob,df,t,U] = clustnum(crds,grps)
  doprob = 0;
  if (nargout>2)
    doprob = 1;
  end;

  [N,p] = size(crds);
  [ng,k] = size(grps);

  if (ng~=N)
    error('  CLUSTNUM: input matrices not compatible');
  end;

  K = zeros(1,k);                       % Numbers of clusters
  for ik = 1:k
    K(ik) = length(unique(grps(:,ik)));
  end;

  [connect,dist] = gabriel(crds,1);     % Find Gabriel graph
  [con,r,c] = trilow(connect);          % List connections
  condist = trilow(dist);               % List corresponding distances

  i = find(con>0);
  r = r(i);
  c = c(i);
  ncon = length(i);

  condist = condist(i);                 % Gabriel distances
  rcondist = ranks(condist);            % Ranks of Gabriel distances

  conid = zeros(ncon,k);                % Connection id: 0=within-cluster, 1=between-cluster

  for ic = 1:ncon                       % For each connection,
    i = r(ic);                            % Get observations
    j = c(ic);
    ig = find(grps(i,:)~=grps(j,:));      % Identify between-cluster dists
    if (length(ig))
      conid(ic,ig) = ones(1,length(ig));
    end;
  end;

  mw = NaN*ones(1,k);
  ma = NaN*ones(1,k);
  if (doprob)
    t =  NaN*ones(1,k);
    U =  NaN*ones(1,k);
    df = NaN*ones(1,k);
    tprob = NaN*ones(1,k);
    Uprob = NaN*ones(1,k);
  end;

  for ik = 1:k                                % For each grp-id vector,
    c = conid(:,ik);                          %   identify within- vs among-cluster ids
    na = sum(c);                                % Number of among-cluster dists
    nw = length(c)-na;                          % Number of within-cluster dists
    
    if (na>0 & nw>0)                            % If there are any clusters,
      mw(ik) = mean(condist(c==0));             %   within-cluster mean dist
      ma(ik) = mean(condist(c==1));             %   among-cluster mean dist
    end;
    if (na>1 & nw>1 & doprob)                          % If there are at least 2 dists within and 2 among,
      [t(ik),tprob(ik),df(ik)] = ttest(condist,c,1);   %   rt-tailed t prob
      [U(ik),Uprob(ik)] =        ttest(rcondist,c,1);  %   rt-tailed t prob of ranks
      
%       tprob(ik) = 6.05e-3 - 1.82e-3*sqrt(N) - 1.34e-3*K(ik) + 0.960*tprob(ik);  % Bias correction
    end;
  end;
  adr = log(ma./mw);

  return;