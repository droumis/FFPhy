% OS_RWNorm: Normalizes the rows of the covariance matrix of normal order statistics  
%            so that the sum of the row elements equals one.  Based on AS 128.4,
%            Davis & Stephens (1978), translated from Fortran.  Called by OrderStats().
%

function v = os_rwnorm(v,n,id)
  small = 1.0e-12;
  nhalf1 = floor((n+1)/2);
  ni = n-1;
  
  for i=2:nhalf1
    s = sum(v(i,i:ni));
    if (id)
      s = s-v(i,i);
    end;
         
    if (abs(s)>small) 
      k=i-1;
      if (id)
        k = i;
      end;
      term = sum(v(i,1:k)) + sum(v(i,(ni+1):n));
      cnst=(1-term)/s;
      m=i;
      if (id)
        m=i+1;
      end;
      nj=n-m+1;
      for j=m:ni
        v(i,j)=v(i,j)*cnst;
        v(j,i)=v(i,j);
        v(ni,nj)=v(i,j);
        v(nj,ni)=v(i,j);
        nj=nj-1;
      end;
    end;
    ni=ni-1;
  end;

  return;

