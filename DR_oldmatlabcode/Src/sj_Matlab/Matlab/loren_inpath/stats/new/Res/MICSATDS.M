% MICSATDS: Calculate genetic distances based on microsatellite loci.
%
%     Usage: dist = micsatds(grps,nloci,freq)
%
%           grps =  [n x 1] vector of group-membership identifiers.
%           nloci = scalar indicating number of loci.
%           freq =  [totallele x 5] matrix of allele frequencies:
%                     col 1 = group identifier
%                         2 = locus identifier
%                         3 = allele identifier (repeat score)
%                         4 = allele frequency within group
%                         5 = standard error of frequency
%           -----------------------------------------------------------------
%           dist =  [k x k] symmetric matrix of intergroup genetic distances.
%

% RE Strauss, 3/11/97

function dist = micsatds(grps,nloci,freq)
  grpid = uniquef(grps);               % Group identifiers
  ngrps = length(grpid);              % Number of groups

  dist = zeros(ngrps,ngrps);

  for k1 = 1:(ngrps-1)                % Get data for group k1
    indx = find(freq(:,1)==grpid(k1));
    Xk1 = freq(indx,:);

    for k2 = (k1+1):ngrps             % Get data for group k2
      indx = find(freq(:,1)==grpid(k2));
      Xk2 = freq(indx,:);

      for c = 1:nloci                 % For each locus,
        indx = find(Xk1(:,2)==c);       % Stash locus data for each grp
        c1 = Xk1(indx,:);
        n1 = size(c1,1);                % Numbers of alleles
        indx = find(Xk2(:,2)==c);
        c2 = Xk2(indx,:);
        n2 = size(c2,1);

        delta = 0;
        for a1 = 1:n1                   % All grp1 alleles vs all grp2 alleles
          f1 = c1(a1,4);                  % Allele freq
          for a2 = 1:n2
            f2 = c2(a2,4);
            delta = delta + f1*f2*(c1(a1,3)-c2(a2,3))^2;
          end;  % for a2=1:n2
        end;  % For a1=1:n1
        dist(k1,k2) = dist(k1,k2) + delta/nloci;
      end;  % For c=1:nloci
      dist(k2,k1) = dist(k1,k2);

    end;  % For k2=2:ngrps
  end;  % For k1=1:(ngrps-1)

  return;
