% BETACF: Evaluates the continued fraction for the incomplete beta function.
%
%     Usage: h = betacf(a,b,x)
%
%        for real values a,b,x, and h

%     Based on the modified Lentz's algorithm of Press et al (1992),
%     "Numerical Recipes in C" (2nd edition), section 6.4.

% RE Strauss, 6/12/93

function h = betacf(a,b,x)
   MAXIT = 100;                     % Maximum iterations
   EPS = 3e-7;
   FPMIN = 1e-30;

   qab = a+b;
   qap = a+1;
   qam = a-1;
   c = 1;
   d = 1 - (qab*x/qap);
   if (abs(d) < FPMIN)
      d = FPMIN;
   end;
   d = 1/d;
   h = d;

   for m=1:MAXIT
      m2 = 2*m;
      aa = m*(b-m)*x/((qam+m2)*(a+m2));
      d = 1+(aa*d);
      if (abs(d) < FPMIN)
         d = FPMIN;
      end;
      c = 1+(aa/c);
      if (abs(c) < FPMIN)
         c = FPMIN;
      end;
      d = 1/d;
      h = h*d*c;
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
      d = 1+(aa*d);
      if (abs(d) < FPMIN)
         d = FPMIN;
      end;
      c = 1+(aa/c);
      if (abs(c) < FPMIN)
         c = FPMIN;
      end;
      d = 1/d;
      del = d*c;
      h = h*del;
      if (abs(del-1) < EPS)
         break;
      end;
   end;

   if (m == MAXIT)
      error('   a or b too bit, or MAXIT too small, in betacf');
   end;
   return;
