% Rankits:  Estimates expected normal order statistics, the expected value of the rth
%           largest value in a sample of size n from a normal distribution.
%           Based on Royston (1982), Applied Statistics algorithm AS 177.
%           Estimates are accurate to at least 0.0001, and usually to 5 or 6 
%           decimal places.  Results have been validated for up to n=2000.
%
%     Usage: z = rankits(n)
%
%         n =    sample size drawn from a normal distribution.
%         ----------------------------------------------------
%         z =    corresponding column vector of expected values.
%

% Royston, J.P. 1985. Algorithm AS177: expected normal order statistics (exact and 
%   approximate).  Applied Statistics 28:161-165.

% RE Strauss, 10/17/02

function z = rankits(n)
  if (~isintegr(n) | n<=0)
    error('  RANKITS: input value must be a positive integer.');
  end;

  n2 = floor(n/2);
  s = zeros(1,n2);

  epsilon = [0.419885,0.450536,0.456936,0.468488];
  delta1 =  [0.112063,0.121770,0.239299,0.215159];
  delta2 =  [0.080122,0.111348,-0.211867,-0.115049];
  gamma =   [0.474798,0.469051,0.208597,0.259784];
  lambda =  [0.282765,0.304856,0.407708,0.414093];
  bb = -0.283833;
  d =  -0.106136;
  s1 =  0.5641896;
  c1 = [9.5,28.7,1.9,0,-7.0,-6.2,-1.6];
  c2 = [-6195,-9569,-6728,-17614,-8278,-3570,1075];
  c3 = [93380,175160,410400,2157000,2376000,2065000,2065000];
  
  if (n==1)
    s = 0;
  elseif (n==2)
    s = s1;
  else
    k = min([n2,3]);
    for i = 1:k
      Q = (i-epsilon(i))/(n+gamma(i));
      Qlambda = Q.^lambda(i);
      correc = 0;
      if (i*n==4)
        correc = 1.9e-5;
      elseif (i<=7 & n<=20)
        an = 1/(n*n);
        correc = (c1(i)+an*(c2(i)+an*c3(i)))*1.0e-6;
      end;
      s(i) = Q+Qlambda*(delta1(i)+Qlambda*delta2(i))/n - correc;
    end;
    if (n2~=k)
      for i = 4:n2
        lambda1 = lambda(4)+bb/(i+d);
        Q = (i-epsilon(4))/(n+gamma(4));
        Qlambda = Q.^lambda1;
        correc = 0;
      if (i*n==4)
        correc = 1.9e-5;
      elseif (i<=7 & n<=20)
        an = 1/(n*n);
        correc = (c1(i)+an*(c2(i)+an*c3(i)))*1.0e-6;
      end;
       s(i) = Q+Qlambda*(delta1(4)+Qlambda*delta2(4))/n - correc;
      end;
    end;
    s = norminv(s);
  end;
  
  z = zeros(n,1);
  r = 1:length(s);
  z(r) = -abs(s);
  z(n-r+1) = -z(r);
  
  return;
  