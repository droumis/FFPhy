% OrderStats: Estimates the means (rankits) and the correlation-covariance matrix of 
%             normal order statistics (ranked values) for n observations, using the
%             methods of Royston (1985) and of Davis & Stephens (1978; AS 128).
%
%     Usage: [m,V] = orderstats(n)
%
%         n = number of observations drawn from a standard-normal distribution.
%         -------------------------------------------------------------------------------
%         m = [1 x n] vector of expected values of the ranked observations.
%         V = [n x n] square symmetric matrix of corresponding variances and covariances.
%

% Royston, J.P. 1985. Algorithm AS177: expected normal order statistics (exact and 
%   approximate).  Applied Statistics 28:161-165.
% Davis, CS & MA Stephens. 1978. Algorithm AS128: approximating the covariance matrix
%   of normal order statistics.  Applied Statistics 27:206-212.

% RE Strauss, 12/9/02

function [m,v] = orderstats(n)
  m = rankits(n);
  ex1 = m(n);
  ex2 = m(n-1);
  sm2 = m'*m;
  
  V = zeros(n,n);

  n2=n+2;
  n22=n2*n2;
  n23=n22*n2;
  nh=floor((n+1)/2);

  ni=n;
  for i=1:nh
    pr=i/(n+1);
    qr=1-pr;
    xr = normcdf(pr);
    [dxr,d2xr,d3xr,d4xr,d5xr] = os_der(xr);
    for j=i:ni
      if (i~=j)
        ps=j/(n+1);
        xs=normcdf(ps);
        [dxs,d2xs,d3xs,d4xs,d5xs] = os_der(xs);
        v(i,j) = os_cov(dxr,d2xr,d3xr,d4xr,d5xr,pr,qr,...
                    dxs,d2xs,d3xs,d4xs,d5xs,ps,n2,n22,n23);
        v(j,i)=v(i,j);
      else
        v(i,j) = os_var(dxr,d2xr,d3xr,d4xr,d5xr,pr,qr,n2,n22,n23);
      end;
    end;
    ni=ni-1;
  end;

  nj = n;
  for i = 2:n
    njm1 = nj - 1;
    im1 = i - 1;
    for j = nj:n
      v(i,j) = v(im1, njm1);
      im1 = im1 - 1;
    end;
    nj = nj - 1;
  end;

  v(1,1) = os_v11(n);
  v(n,n)=v(1,1);
  v(1,2)=v(1,1)+ex1*(ex1-ex2)-1;
  v(2,1)=v(1,2);
  v(n,n-1)=v(1,2);
  v(n-1,n)=v(1,2);
  
  if (n==2)
    retun;
  end;
    
  s = sum(v(1,3:n));  
  cnst=(1-v(1,1)-v(1,2))/s;
  nj=n-2;
  for j=3:n
    v(1,j)=v(1,j)*cnst;
    v(j,1)=v(1,j);
    v(n,nj)=v(1,j);
    v(nj,n)=v(1,j);
    nj=nj-1;
  end;

  v = os_rwnorm(v,n,0);
      
  s=0;
  for k=1:n
    if (k~=2 & k~=n-1)
      s=s+v(k,k);
    end;
  end;
  v(2,2)=0.5*(n-sm2-s);
  v(n-1,n-1)=v(2,2);

  v = os_rwnorm(v,n,1);

  return;
  