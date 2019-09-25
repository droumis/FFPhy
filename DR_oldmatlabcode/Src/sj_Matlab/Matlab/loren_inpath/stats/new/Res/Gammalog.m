% GAMMALOG: Logarithm of the gamma function.  Uses Lanczos-type approximation 
%           to ln(gamma) for z > 0.  Accurate to about 14 significant digits 
%           except for small regions in the vicinity of 1 and 2.
%
%     Usage: gl = gammalog(z)
%

% Lanczos, C. 1964. A precision approximation of the gamma function.
%   J. SIAM Numer. Anal. B1:86-96.

% RE Strauss, 3/2/01

function gl = gammalog(z)
  a = [ 0.9999999999995183; 676.5203681218835;    % Constants
       -1259.139216722289;  771.3234287757674,
       -176.6150291498386;  12.50734324009056,
       -0.1385710331296526; 0.9934937113930748e-05,
        0.1659470187408462e-06];
  lnsqrt2pi = log(sqrt(2*pi));

  [r,c] = size(z);
  gl = zeros(size(z));

  for i = 1:r
    for j = 1:c
      zz = z(i,j);
      g = 0;
      t = zz + 7;
      for k = 9:-1:2
        g = g + a(k)/t;
        t = t-1;
      end;
      g = g + a(1);
      gl(i,j) = log(g) + lnsqrt2pi - (zz+6.5) + (zz-0.5)*log(zz+6.5);
    end;
  end;

  return;