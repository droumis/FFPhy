% WRIGHTSF: Based on the results of fitting a single-factor general model using Wright's 
%           method, investigates the structure of the residual covariances to find 
%           plausible secondary factors, based on UPGMA cluster analysis of pairwise 
%           residual covariances.
%
%     Usage: [gen_factor,gen_r2a,sec_vars,prime_factors,sec_factors,sec_r2a] = ...
%                 wrightsf(covmat,{maxsecvar},{allsecfact})
%
%           covmat =        [p x p] covariance or correlation matrix.
%           maxsecvar =     optional maximum number of variables to be included in 
%                             secondary factor [default = floor(p/2)].
%           allsecfact =    optional boolean flag indicating whether (=1) or not (=0) all 
%                             secondary factors examined are to be return; default (=0) 
%                             is that only secondary factors that increase the adjusted 
%                             R^2 are returned.
%           ----------------------------------------------------------------------------
%           gen_factor =    [p x 1] vector of general single-factor coefficients.
%           gen_r2a =       adjusted R^2 (goodness of fit) for general factor.
%           sec_vars =      [s x maxsecvar] matrix of combinations of variables included 
%                             in the 's' secondary factors examined, sequenced by 
%                             decreasing adjusted R^2.
%           prime_factors = [s x p] matrix of corresponding primary-factor coefficients.
%           sec_factors =   [s x maxsecvar] matrix of corresponding secondary-factor 
%                             coefficients.
%           sec_r2a =       [s x 1] vector of corresponding adjusted R^2 values for full 
%                             models (primary + secondary factors).
%           

function [gen_factor,gen_r2a,sec_vars,prime_factors,sec_factors,sec_r2a] = ...
            wrightsf(covmat,maxsecvar,allsecfact)

  if (nargin < 2)
    maxsecvar = [];
  end;
  if (nargin < 3)
    allsecfact = [];
  end;

  nvars = size(covmat,1);                       % Number of variables
  ncovars = 0.5*nvars*(nvars-1);                % Number of covariances

  if (isempty(maxsecvar))                       % Default input arguments
    maxsecvar = floor(nvars/2);
  end;
  if (isempty(allsecfact))
    allsecfact = 0;
  end;

  [gen_factor,gen_r2a,resid_sf] = wright(covmat);  % Fit single general-factor model

  [rc,row,col] = trilow(resid_sf);              % Sort abs(residual covars) in decreasing sequence
  [rcx,rc,row,col] = sortmat(-abs(rc),rc,row,col);

  sec_vars = zeros(ncovars,maxsecvar);          % Allocate output matrices
  prime_factors = zeros(ncovars,nvars);
  sec_factors = zeros(ncovars,maxsecvar);
  sec_r2a  = zeros(ncovars,1);

  i = 0;
  while (i < ncovars)                           % Fit all pairwise combinations of variables
    i = i+1;                                    %   as secondary factors
    v1 = row(i);
    v2 = col(i);

    secnd = zeros(nvars,1);
    secnd([v1 v2]) = [1; 1];
    sec_vars(i,1:2) = [v1 v2];

    [f,r2] = wright(covmat,secnd);

    prime_factors(i,:) = f(:,1)';
    s = f(:,2)';
    sec_factors(i,1:2) = s(abs(s)>0);
    sec_r2a(i) = abs(r2(2));
  end;

  dist = zeros(nvars,nvars);            % Create distance matrix from pairwise 1-R^2 values
  row = sec_vars(:,1);
  col = sec_vars(:,2);
  for i = 1:ncovars
    d = 1-sec_r2a(i);
    dist(row(i),col(i)) = d;
    dist(col(i),row(i)) = d;
  end;

  [topology,support] = upgma(dist);     % Cluster analysis of distance matrix

  sn = sum((support>0)')';              % Numbers of vars in each cluster
  i = find(sn<3 | sn>maxsecvar);        % Delete clusters (rows) having too few or too many
  support(i,:) = [];                    %   vars
  [nclust,c] = size(support);
  if (c > maxsecvar)                    % Reduce to first 'maxsecvar' columns
    support(:,(maxsecvar+1):c) = [];
  end;

  r2_clust = zeros(nclust,1);           % Allocate matrices for cluster-specified
  pf = zeros(nclust,nvars);             %   secondary factors
  sf = zeros(nclust,maxsecvar);

  for i = 1:nclust                      % Cycle thru clusters, specifying vars as
    v = support(i,:);                   %   secondary factors
    v(v==0) = [];

    secnd = zeros(nvars,1);
    secnd(v) = ones(1,length(v));

    [f,r2] = wright(covmat,secnd);

    pf(i,:) = f(:,1)';
    s = f(:,2)';
    sf(i,1:length(v)) = s(abs(s)>0);
    r2_clust(i) = abs(r2(2));
  end;

  sec_vars = [sec_vars; support];       % Append cluster results to output matrices
  prime_factors = [prime_factors; pf];
  sec_factors = [sec_factors; sf];
  sec_r2a = [sec_r2a; r2_clust];

  s = sum(support);                     % Remove columns of secondary-factor matrices
  if (any(s==0))                        %   having no values
    i = find(s==0);
    sec_vars(:,i) = [];
    sec_factors(:,i) = [];
  end;

  [sec_r2a,sec_vars,prime_factors,sec_factors] = ...    % Sort matrices by decreasing R^2
      sortmat(-sec_r2a,sec_vars,prime_factors,sec_factors);
  sec_r2a = -sec_r2a;

  if (~allsecfact)                      % Reduce secondary factors to those that
    i = (sec_r2a >= gen_r2a);           %   increase the R^2 over the general factor
    sec_vars = sec_vars(i,:);
    prime_factors = prime_factors(i,:);
    sec_factors = sec_factors(i,:);
    sec_r2a = sec_r2a(i);
  end;

  return;

