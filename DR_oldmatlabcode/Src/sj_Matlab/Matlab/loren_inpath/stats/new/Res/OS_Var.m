% OS_Var: Calculates David-Johnson approximation for the variance of the rth 
%         largest order statistic from the normal distribution for sample size n.  
%         Based on AS 128.2, Davis & Stephens (1978), translated from Fortran.  
%         Called by OrderStats().
%

function v = os_var(dxr,d2xr,d3xr,d4xr,d5xr,pr,qr,rn2,rn22,rn23)

  dxr2=dxr*dxr;
  prqr=pr*qr;
  v=prqr*dxr2/rn2;
  qrmpr=qr-pr;
  d2xr2=d2xr*d2xr;
  v = v+prqr/rn22*(2*qrmpr*dxr*d2xr+prqr*(dxr*d3xr+0.5*d2xr2));
  v = v+prqr/rn23*(-2*qrmpr*dxr*d2xr...
    +(qrmpr*qrmpr-prqr)*(2*dxr*d3xr+1.5*d2xr2)...
    +prqr*qrmpr*((5/3)*dxr*d4xr+3*d2xr*d3xr)...
    +0.25*prqr*prqr*(dxr*d5xr+2*d2xr*d4xr+(5/3)*d3xr*d3xr));
    
  return;
  
