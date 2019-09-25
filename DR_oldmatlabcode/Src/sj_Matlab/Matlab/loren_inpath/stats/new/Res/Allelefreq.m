% ALLELEFREQ: Calculates allele frequencies.  Ignores missing data.
%
%     Usage: [freq,avgfreq] = allelefreq(alleles,popl)
%
%           alleles = [n x 2c] matrix of allele identifiers for n obs and c loci.
%           popl =    [n x 1] vector of group-membership identifiers for k groups.
%           ----------------------------------------------------------------------
%           freq =    [totallele x 6] matrix of allele frequencies
%                       col 1 = group identifier
%                           2 = locus identifier
%                           3 = allele identifier
%                           4 = allele count
%                           5 = within-group allele frequency 
%                           6 = standard error of frequency
%           avgfreq = [totallele x 4] matrix of frequencies pooled across groups
%                       col 1 = locus identifier
%                           2 = allele identifier
%                           3 = mean allele frequency                       
%                           4 = standard error of mean frequency
%

% RE Strauss, 3/11/97
%   9/20/99 - update handling of null input arguments.
%   5/13/00 - added standard errors of frequencies and heterozygosity estimates; 
%             renamed from micsatfr.
%   5/16/00 - added expected heterozygosity estimates; added allele counts to 
%               frequency matrix.
%   5/25/00 - separated heterozygosity estimates into function heterozyg().

function [freq,avgfreq] = allelefreq(alleles,popl)
  popl = popl(:);                     % popl should be a column vector
  grpid = uniquef(popl);              % Group identifiers
  npopl = length(grpid);              % Number of groups (k)
  [nobs,i] = size(alleles);           % Number of observations and loci (c)
  nloci = floor(i/2);

  if (~isvector(popl))
    error('  ALLELEFREQ: grouping vector must be a vector.');
  end;
  if (length(popl) ~= nobs)
    error('  ALLELEFREQ: allele matrix and grouping vector not compatible.');
  end;

  get_avgfreq = 0;

  if (nargout >= 2)                     % Set flags for output operations
    get_avgfreq = 1;
  end;

  % Allele frequencies
  %   freq = [totallele x 6] matrix of allele frequencies
  %            col 1 = group identifier
  %                2 = locus identifier
  %                3 = allele identifier
  %                4 = allele count
  %                5 = within-group allele frequency 
  %                6 = standard error of frequency

  freq = [];

  for k = 1:npopl 
    ig = find(popl==grpid(k));          % Locate data for current group

    for c = 1:nloci
      Xkc = alleles(ig,(2*c-1):(2*c));  % Stash current locus
      Xkc = Xkc(:);                     % Concatenate the two allele vectors
      indx = find(~isfinite(Xkc));      % Delete missing data
      if (~isempty(indx))
        Xkc(indx) = [];
      end;
      nk = length(Xkc);                 % Sample size for this locus

      if (nk>0)
        allid =  uniquef(Xkc);          % Unique allele ids for current locus
        allid = sort(allid);            % Sort into ascending sequence
        nallid = length(allid);         % Number of unique alleles

        f = zeros(nallid,1);            % Corresponding allele frequencies
        sef = zeros(nallid,1);          %   standard errors
        cnt = zeros(nallid,1);          %   counts
        for a = 1:nallid
          indx = find(Xkc==allid(a));
          na = length(indx);
          ff = na / nk;
          f(a) = ff;
          sef(a) = sqrt(ff*(1-ff)/nk);
          cnt(a) = na;
        end;

        freq = [freq; grpid(k)*ones(nallid,1), c*ones(nallid,1), allid, cnt, f, sef];
      end;
    end;  % For loci
  end;  % For npopl


  % Average allele frequencies across all groups
  %   avgfreq = [totallele x 4] matrix of frequencies pooled across groups
  %               col 1 = locus identifier
  %                   2 = allele identifier
  %                   3 = mean allele frequency                       
  %                   4 = standard error of mean frequency

  avgfreq = [];       
                
  if (get_avgfreq)
    for c = 1:nloci
      Xkc = alleles(:,(2*c-1):(2*c));   % Stash current locus
      Xkc = Xkc(:);                     % Concatenate the two allele vectors
      indx = find(~isfinite(Xkc));      % Delete missing data
      if (~isempty(indx))
        Xkc(indx) = [];
      end;
      nk = length(Xkc);                 % Sample size for this locus

      if (nk>0)
        allid =  uniquef(Xkc);          % Unique allele ids for current locus
        allid = sort(allid);            % Sort into ascending sequence
        nallid = length(allid);         % Number of unique alleles

        f = zeros(nallid,1);            % Corresponding allele frequencies
        sef = zeros(nallid,1);          %   and standard errors
        for a = 1:nallid
          indx = find(Xkc==allid(a));
          ff = length(indx) / nk;
          f(a) = ff;
          sef(a) = sqrt(ff*(1-ff)/nk);
        end;
  
        avgfreq = [avgfreq; c*ones(nallid,1), allid, f, sef];
      end;
    end;  % For loci
  end;

  return;
