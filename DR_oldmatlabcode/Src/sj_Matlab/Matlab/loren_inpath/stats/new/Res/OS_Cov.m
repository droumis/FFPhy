% OS_Cov: Calculates David-Johnson approximation for the covariance between the 
%         rth and sth largest order statistics from the normal distribution for  
%         sample size n. Based on AS 128.3, Davis & Stephens (1978), translated 
%         from Fortran. Called by OrderStats().
%

function cv = os_cov(dxr,d2xr,d3xr,d4xr,d5xr,pr,qr,...
                    dxs,d2xs,d3xs,d4xs,d5xs,ps,rn2,rn22,rn23)
  qs=1-ps;
  prqs=pr*qs;
  cv=prqs*dxr*dxs/rn2;
  qrmpr=qr-pr;
  qsmps=qs-ps;
  prqr=pr*qr;
  psqs=ps*qs;
  cv = cv+prqs/rn22*(qrmpr*d2xr*dxs+qsmps*dxr*d2xs+0.5*prqr*d3xr*dxs...
        +0.5*psqs*dxr*d3xs+0.5*prqs*d2xr*d2xs);
  pr2=pr*pr;
  qr2=qr*qr;
  ps2=ps*ps;
  qs2=qs*qs;
  psqr=ps*qr;
  term1=-d2xr*dxs*qrmpr-qsmps*dxr*d2xs+(qrmpr*qrmpr-prqr)*d3xr*dxs;
  term2=(qsmps*qsmps-psqs)*dxr*d3xs+(1.5*qrmpr*qsmps...
         +0.5*psqr-2*prqs)*d2xr*d2xs;
  term3=(5/6)*(prqr*qrmpr*d4xr*dxs+psqs*qsmps*dxr*d4xs)...
         +(prqs*qrmpr+0.5*prqr*qsmps)*d3xr*d2xs;
  term4=(prqs*qsmps+0.5*psqs*qrmpr)*d2xr*d3xs...
         +(1/8)*(pr2*qr2*d5xr*dxs+ps2*qs2*dxr*d5xs);
  term5=0.25*(pr2*qr*qs*d4xr*d2xs+pr*ps*qs2*d2xr*d4xs)+(1/12)...
         *(2*pr2*qs2+3*pr*qr*ps*qs)*d3xr*d3xs;
  cv = cv+prqs/rn23*(term1+term2+term3+term4+term5);
      
  return;
