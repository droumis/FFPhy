% FACTORS: Estimates loadings of a set of secondary factors of a covariance
%          matrix, given a primary factor and a corresponding boolean matrix
%          of specifying variables included in the secondary factor(s).
%          Secondary factors by default are orthogonal (with the primary factor
%          and with one another), but this constraint can be relaxed.
%          Also returns the residual covariance matrix.
%          Error handling is done by function WRIGHT.
%
%     Syntax: [fact,rc] = factors(c,pfactor,secnd,orthog,suppress)
%
%         c =       [n x n] residual covariance or correlation matrix.
%         pfactor = n-length vector of loadings of primary factor.
%         secnd =     [n x s] matrix of boolean values (T/F = 1/0)
%                     indicating by 1's the submatrices of variables
%                     to be included in the secondary factors (s of them).
%         orthog =  flag indicating whether secondary factors are to be:
%                     0 - oblique to primary factor and one another;
%                     1 - orthogonal to primary factor but oblique to other
%                           secondary factors;
%                     2 - orthogonal to primary factor and all other secondary
%                           factors;
%                     [default = 2].
%         suppress = boolean flag indicating whether warning messages
%                     are to be suppressed (invoked for bootstrapping and
%                     randomization).
%         ----------------------------------------------------------------
%         fact =    [n x s+1] matrix of factor loadings (primary + secondary).
%         rc =      [n x n] matrix of residual covariances.
%

function [fact,rc] = factors(c,pfactor,secnd,orthog,suppress)
  TRUE = 1; FALSE = 0;

  [nvar,p] = size(c);
  [n,nsec] = size(secnd);

  % Estimate loadings from covariance matrix.  Initial guess is that
  % loadings are equal to the root-mean-square covariance
  % (including variance), by column, orthogonal to other factors.

  % First pass: optimize secondary factors separately with respect to
  % primary factor only.

  fact = zeros(nvar,nsec);
  init_load = [];

  for s = 1:nsec                      % Initial loadings for secondary factors
    secvect = find(secnd(:,s));           % Transform current factor into var list

    if (length(secvect)==2)             % 2 variables: exact solution
      cc = c(secvect(1),secvect(2));      % Covar between the 2 variables

      if (orthog==0)                      % If factor can be oblique to primary,
        f = ones(2,1) * sqrt(abs(cc));      % Each loading is sqrt of covar
        if (cc<0)                           % If covar is negative, one loading
          f(2) = -f(2);                     %   is negative
        end;
      else                                % If should be orthogonal to primary,
        a1 = pfactor(secvect(1));           % Exact solution possible only if:
        a2 = pfactor(secvect(2));           %   primary factor is general,
        caa = -cc*a2/a1;                    %   secondary factor is bipolar.
        if (caa > 0)
          f = ones(2,1);
          f(1) = sqrt(caa);
          f(2) = cc / f(1);
        else                                % But not possible if covar is pos,
          f = ones(2,1) * sqrt(abs(cc));      % so each loading is sqrt of covar
          if (~suppress)
            disp('    Warning: 2-variable secondary factor must be oblique');
            disp('             to primary factor when residual covariance');
            disp('             is positive.');
            suppress = TRUE;
          end;
        end;
      end;

    else                                % >2 variables: estimated solution
      redc = c(secvect,secvect);          % Reduced covar matrix for this factor
      g = min(sqrt(abs(mean(redc)))',sqrt(diag(redc)) );         % Initial loadings
      f = fmins('factorsf',g,[],[],c,pfactor,secnd(:,s),orthog); % Optimize this factor
    end;

    fact(secvect,s) = f;                % Insert into factor-loadings matrix
    init_load = [init_load;f];          % Append onto initial-loadings vector
  end;

  % Optional second pass: simultaneously optimize 2+ secondary factors if they
  % are intercorrelated.  Secondary factors of two variables having loadings
  % of the same sign are omitted from second-pass optimization because they
  % must be oblique to any primary factor.

  if (orthog)
    ns = sum(secnd);                      % Omit sec factors with 2 vars of same sign

    omit = [];
    keep = [];
    for s = 1:nsec
      if (ns(s) > 2)
        keep = [keep s];
      else
        if (abs(sum(sign(fact(:,s)))))
          omit = [omit s];
        else
          keep = [keep s];
        end;
      end;
    end;

    if (nsec-length(omit) > 1)          % Simultaneously optimize >1 factor
      if (length(omit)>0)                 % If 2-var factors to be omitted,
        nsec_copy = nsec;                 %   Make copies of full matrices
        sec_copy = secnd;
        fact_copy = fact;

        nsec = nsec - length(omit);       %   Delete oblique 2-var factors
        secnd(:,omit) = [];
        fact(:,omit) = [];

        factors_omitted = TRUE;
      else
        factors_omitted = FALSE;
      end;

      vc = abs(vectcorr(fact));           % Vector corrs among secondary factors
      sumcorr = sum(trilow(vc));

      if (sumcorr > eps)                  % Optimize if any corrs > 0
        f = fmins('factorsf',init_load,[],[],c,pfactor,secnd,orthog);
      end;

      low = 1;                            % Unpack loadings vector into matrix
      for s = 1:nsec                      %   of secondary-factor loadings
        secvect = find(secnd(:,s));
        high = low+length(secvect)-1;
        fact(secvect,s) = f(low:high);
        low = high+1;
      end;
    end;

    if (factors_omitted)                % If 2-var factors have been omitted,
      nsec = nsec_copy;                 %   restore them
      secnd = sec_copy;
  
      fact_copy(:,keep) = fact;
      fact = fact_copy;
    end;
  end;  % if orthog

%  r2 = vectcorr(pfactor,fact).^2;     % Independent variances of secondary factors
%  rc = c;                             % Residual covariances
  for s = 1:nsec
    f = fact(:,s);
    ff = f*f';                          % Reproduced covariances
%    for i = 1:nvar
%      if (abs(f(i))>eps)
%       rc(i,i) = rc(i,i) - ff(i,i);    % Fully reproduce diagonal elements
%        ff(i,i) = 0;                    % Conditionally reproduce off-diag elements
%      end;
%    end;
%    rc = rc - ff*(1-r2(s));
    rc = c - ff;
  end;

  fact = [pfactor fact];              % Return primary + secondary factor loadings

  return;
