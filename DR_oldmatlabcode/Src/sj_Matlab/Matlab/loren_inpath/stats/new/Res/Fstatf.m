% FSTATF: Objective function for FSTAT: estimation of Wright's F-statistics 
%         (FST [theta], FIS [f], FIT [F]).
%
%
%     Usage: fvect = fstatf(alleles,popl,nu1,nu2)
%
%         alleles =   [N x 2*c] matrix of allele identities for N specimens 
%                       (grouped into r populations) and c loci.
%                       Alleles assumed to be coded 1,2,...,m
%                       Missing values should be coded as zeros.
%         popl =      grouping vector (length N) containing ids of populations 
%                       for all specimens.
%         nu1, nu2 =  unused arguments passed by bootstrp().
%         ---------------------------------------------------------------------
%         fvect =     row vector (length (c+1)*3) of [Fis Fit Fst] values.
%

% Weir, BS and CC Cockerham. 1984. Estimating F-statistics for the analysis of 
%   population structure.  Evolution 38:1358-1370.

% RE Strauss, 8/22/00

function fvect = fstatf(alleles,popl,nu1,nu2)
  nloci = size(alleles,2)/2;
  
  locus = [];
  aa = [];
  ab = [];
  ac = [];

  for loc = 1:nloci                     % For all loci,
    a1 = alleles(:,2*loc-1);              % Pull out alleles
    a2 = alleles(:,2*loc);

    i = find(~isfinite(a1) | ~isfinite(a2)); % Delete specs with missing data
    a1(i) = [];
    a2(i) = [];
    pop = popl;
    pop(i) = [];

    if (~isempty(popl))
      [ua,fa] = uniquef([a1;a2]);           % Unique alleles and freqs
      m = length(ua);                       % Number of alleles
      [popid,n] = uniquef(pop);             % Unique popl identifiers & sizes
      r = length(popid);                    % Number of populations

      if (m>1 & r>1)                        % If locus polymorphic & have >1 population,
        rm1 = r-1;
        nbar = mean(n);                       % Mean popl size
        ntot = sum(n);

        freq = zeros(r,1);                    % Alloc allele freqs per popl
        het = zeros(r,1);                     % Alloc heterozygote freqs per popl

        for all = 1:m                         % For all alleles,
          aid = ua(all);                        % Current allele identifier
          count = fa(all);                      % Total allele count
          meanfreq = count/sum(fa);             % Weighted mean allele freq across popls

          for p = 1:r                           % For all populations,
            i = find(pop==popid(p));              % Find specimens
            a1i = a1(i);
            a2i = a2(i);

            freq(p) = sum([a1i;a2i]==aid)/(2*n(p));     % Allele freq

            j = find(a1i==aid | a2i==aid);
            if (~isempty(j))
              a1ij = a1i(j);
              a2ij = a2i(j);
              het(p) = sum((a1ij==aid & a2ij~=aid) | (a1ij~=aid & a2ij==aid))/n(p);
            else
              het(p) = 0;
            end;
          end;

          s2a = sum(n.*((freq-meanfreq).^2))/(rm1*nbar);
          Ha =  sum(n.*het)/ntot;
          nc =  (ntot-(sum(n.*n)/ntot))/rm1;
          pq =  meanfreq*(1-meanfreq);
          nbarm1 = nbar-1;
          nbarmnc = nbar-nc;

          a = (s2a - (pq - rm1*s2a/r - Ha/4)/nbarm1) * nbar/nc;
          b = (pq - rm1*s2a/r - Ha*(2*nbar-1)/(4*nbar)) * nbar/nbarm1;
          c = Ha/2;

          if (abs(b) > eps)
            locus = [locus; loc];
            aa = [aa; a];
            ab = [ab; b];
            ac = [ac; c];
          end;
        end;
      end;
    end;
  end;

  fval = NaN*ones(nloci+1,3);
  for loc = 1:nloci
    i = find(locus==loc);
    if (~isempty(i))
      sumabc = sum(aa(i))+sum(ab(i))+sum(ac(i));
      lfit = 1 - sum(ac(i))/sumabc;
      lfst = sum(aa(i))/sumabc;
      lfis = 1 - sum(ac(i))/(sum(ab(i))+sum(ac(i)));
      fval(loc,:) = [lfis lfit lfst];
    end;
  end;
 
  sumabc = sum(aa)+sum(ab)+sum(ac);
  fit = 1 - sum(ac)/sumabc;
  fst = sum(aa)/sumabc;
  fis = 1 - sum(ac)/(sum(ab)+sum(ac));
  fval(nloci+1,:) = [fis fit fst];

  fvect = fval(:)';

  return;
