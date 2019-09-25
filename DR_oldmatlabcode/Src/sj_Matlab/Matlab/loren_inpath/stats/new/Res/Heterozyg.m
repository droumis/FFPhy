% HETEROZYG: Calculates observed and expected heterozygosity, and their unbiased
%            standand errors, from data matrix.  Ignores missing data.
%            (Source: Weir, BS. 1996. Genetic Data Analysis II. Sinauer)
%
%     Usage: [hetzyg,avghetzyg,exphetzyg] = heterozyg(alleles,popl)
%
%           alleles = [n x 2c] matrix of allele identifiers for n obs and c loci.
%           popl =    [n x 1] vector of group-membership identifiers for k groups.
%           ----------------------------------------------------------------------
%           hetzyg =  [k*(c+1) x 4] matrix of within-group heterozygosity estimates
%                       col 1 = group identifier
%                           2 = locus identifier (0 = mean across loci)
%                           3 = within-group heterozygosity
%                           4 = standard error of heterozygosity
%           avghetzyg = [(c+1) x 3] matrix of mean heterozygosity across groups
%                       col 1 = locus identifier (0 = mean across loci)
%                           2 = mean heterozygosity
%                           3 = standard error of mean heterozygosity
%           exphetzyg = [k*(c+1) x 4] matrix of expected heterozygosity estimates
%                       col 1 = group identifier
%                           2 = locus identifier (0 = mean across loci)
%                           3 = expected within-group heterozygosity
%                           4 = standard error of expected heterozygosity
%

% RE Strauss, 3/11/97
%   9/20/99 - update handling of null input arguments.
%   5/13/00 - added standard errors of frequencies and heterozygosity estimates; 
%             renamed from micsatfr.
%   5/16/00 - added expected heterozygosity estimates; added allele counts to 
%               frequency matrix.
%   5/25/00 - separated from function allelefreqs.

function [hetzyg,avghetzyg,exphetzyg] = heterozyg(alleles,popl)
  popl = popl(:);                     % popl should be a column vector
  grpid = uniquef(popl);              % Group identifiers
  npopl = length(grpid);              % Number of groups (k)
  [nobs,i] = size(alleles);           % Number of observations and loci (c)
  nloci = floor(i/2);

  if (~isvector(popl))
    error('  HETEROZYG: grouping vector must be a vector.');
  end;
  if (length(popl) ~= nobs)
    error('  HETEROZYG: allele matrix and grouping vector not compatible.');
  end;

  get_exphetzyg = 0;
  get_avghetzyg = 0;

  if (nargout >= 2)                     % Set flags for output operations
    get_avghetzyg = 1;
  end;
  if (nargout >= 3)
    get_exphetzyg = 1;    
  end;


  % Within-group heterozygosity estimates, by locus
  %   hetzyg =  [k*(c+1) x 4] matrix of within-group heterozygosity estimates
  %               col 1 = group identifier
  %                   2 = locus identifier (0 = mean across loci)
  %                   3 = within-group heterozygosity
  %                   4 = standard error of heterozygosity

  hetzyg = zeros(npopl*(nloci+1),4);
  ih = 0;
  hv = zeros(nloci,1);
  hn = zeros(nloci,1);

  for k = 1:npopl                     % For all groups
    ig = find(popl==grpid(k));          % Locate data for current group
    Xk = alleles(ig,:);                 % Isolate current group

    for c = 1:nloci                     % Single-locus heterozygosities
      Xkc = Xk(:,(2*c-1):(2*c));          % Stash current locus

      i = find(~isfinite(sum(Xkc')));     % Delete missing data
      if (~isempty(i))
        Xkc(i,:) = [];
      end;

      n = size(Xkc,1);                    % Number of individuals
      nhet = length(find(Xkc(:,1) ~= Xkc(:,2)));  % Number of heterozygotes
      h = nhet/n;
      varh = h*(1-h)/n;

      hv(c) = h;
      hn(c) = n;

      ih = ih+1;
      hetzyg(ih,:) = [grpid(k) c h sqrt(varh)];
    end;  % For loci

    avgh = (hv'*hn)/sum(hn);              % Weighted mean
    varh = sum(hv.*(1-hv)./hn);           % Single-het component of variance

    for c1 = 1:nloci                      % Double-het component of variance
      for c2 = 1:nloci                    %   at pairs of loci
        if (c2~=c1)
          Xkc = Xk(:,[(2*c1-1) (2*c1) (2*c2-1) (2*c2)]);

          i = find(sum(Xkc'));              % Delete missing data
          if (~isfinite(i))
            Xkc(i,:) = [];
          end;
          nn = size(Xkc,1);                 % Number of individuals for both loci

          if (nn > 0)
            h1 = length(find(Xkc(:,1)~=Xkc(:,2)))/nn;  % Single & double hets
            h2 = length(find(Xkc(:,3)~=Xkc(:,4)))/nn;
            h12 = length(find(Xkc(:,1)~=Xkc(:,2) & Xkc(:,3)~=Xkc(:,4)))/nn;
            varh = varh + (h12 - h1*h2);
          end;
        end;
      end;
    end;

    varh = abs(varh)/(nloci^2);
    ih = ih+1;
    hetzyg(ih,:) = [grpid(k) 0 avgh sqrt(varh)];
  end;  % For groups


% Mean heterozygosity, by locus across groups
%   avghetzyg = [(c+1) x 3] matrix of mean heterozygosity across groups
%               col 1 = locus identifier (0 = mean across loci)
%                   2 = mean heterozygosity
%                   3 = standard error of mean heterozygosity

  if (get_avghetzyg)
    avg_hetzyg = zeros(nloci+1,3);
    ih = 0;
    hv = zeros(nloci,1);
    hn = zeros(nloci,1);

    for c = 1:nloci                     % Single-locus heterozygosities
      Xkc = alleles(:,(2*c-1):(2*c));          % Stash current locus

      i = find(~isfinite(sum(Xkc')));     % Delete missing data
      if (~isempty(i))
        Xkc(i,:) = [];
      end;

      n = size(Xkc,1);                    % Number of individuals
      nhet = length(find(Xkc(:,1) ~= Xkc(:,2)));  % Number of heterozygotes
      h = nhet/n;
      varh = h*(1-h)/n;

      hv(c) = h;
      hn(c) = n;

      ih = ih+1;
      avghetzyg(ih,:) = [c h sqrt(varh)];
    end;  % For loci

    avgh = (hv'*hn)/sum(hn);              % Weighted mean
    varh = sum(hv.*(1-hv)./hn);           % Single-het component of variance

    for c1 = 1:nloci                      % Double-het component of variance
      for c2 = 1:nloci                    %   at pairs of loci
        if (c2~=c1)
          Xkc = alleles(:,[(2*c1-1) (2*c1) (2*c2-1) (2*c2)]);

          i = find(sum(Xkc'));              % Delete missing data
          if (~isfinite(i))
            Xkc(i,:) = [];
          end;
          nn = size(Xkc,1);                 % Number of individuals for both loci

          if (nn > 0)
            h1 = length(find(Xkc(:,1)~=Xkc(:,2)))/nn;  % Single & double hets
            h2 = length(find(Xkc(:,3)~=Xkc(:,4)))/nn;
            h12 = length(find(Xkc(:,1)~=Xkc(:,2) & Xkc(:,3)~=Xkc(:,4)))/nn;
            varh = varh + (h12 - h1*h2);
          end;
        end;
      end;
    end;

    varh = abs(varh)/(nloci^2);
    ih = ih+1;
    avghetzyg(ih,:) = [0 avgh sqrt(varh)];
  end;


  % Expected within-group heterozygosity estimates, by locus
  %   hetzyg =  [k*(c+1) x 4] matrix of expected heterozygosity estimates
  %               col 1 = group identifier
  %                   2 = locus identifier (0 = mean across loci)
  %                   3 = expected heterozygosity
  %                   4 = standard error of heterozygosity

  if (get_exphetzyg)
    freq = allelefreq(alleles,popl);    % Get allele frequencies

    exphetzyg = zeros(npopl*(nloci+1),4);
    ih = 0;
    hv = zeros(nloci,1);
    hn = zeros(nloci,1);

    for k = 1:npopl                     % For all groups
      for c = 1:nloci
        i = find(freq(:,1)==grpid(k) & freq(:,2)==c); % Isolate frequencies
        n = sum(freq(i,4));
        f = freq(i,5);
        h = 1 - f'*f;                     % Heterozygosity under HW frequencies
        varh = h*(1-h)/n;

        hv(c) = h;
        hn(c) = n;

        ih = ih+1;
        exphetzyg(ih,:) = [grpid(k) c h sqrt(varh)];
      end;

      avgh = (hv'*hn)/sum(hn);              % Weighted mean
      varh = sum(hv.*(1-hv)./hn);           % Single-het component of variance
      varh = abs(varh)/(nloci^2);

      ih = ih+1;
      exphetzyg(ih,:) = [grpid(k) 0 avgh sqrt(varh)];
    end;  % For groups
  end;

  return;
