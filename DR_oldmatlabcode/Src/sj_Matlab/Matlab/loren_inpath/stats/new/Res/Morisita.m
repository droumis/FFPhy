% MORISITA: Morisita's measure of dissimilarity between two communities (sites), 
%           based on counts of individuals of composite species.
%           Note: can produce negative distances for very small sample sizes.
%
%     Syntax: dist = morisita(abund)
%
%           abund = [N x S] matrix of absolute abundance (counts of 
%                     individuals) for N species across S sites.
%           ---------------------------------------------------------------------
%           dist =  [S x S] symmetric distance matrix.
%

% RE Strauss, 2/26/97
%   6/29/00 - correct problem with zero lambdas.

function dist = morisita(abund)
  [n,s] = size(abund);

  sim = ones(s,s);
  sumabund = sum(abund);

  for i = 1:(s-1)                   % All possible pairs of sites (columns)
    ni = sumabund(i);
    if (ni < 2)
      sim(i,:) = zeros(1,s);
      sim(:,i) = zeros(s,1);
      sim(i,i) = 1;
    else
      lambdai = sum(abund(:,i).*(abund(:,i)-1))/(ni*(ni-1));
      for j = (i+1):s
        nj = sumabund(j);
        sim(i,j) = 0;
        sim(j,i) = 0;
        if (nj > 1)
          lambdaj = sum(abund(:,j).*(abund(:,j)-1))/(nj*(nj-1));
          crossprod = sum(abund(:,i).*abund(:,j));
          if (lambdai>0 & lambdaj>0)
            sim(i,j) = 2*crossprod/((lambdai+lambdaj)*ni*nj);
            sim(j,i) = sim(i,j);
          end;
        end;
      end;
    end;
  end;

  dist = 1-sim;                    % Convert to distance matrix

  return;
