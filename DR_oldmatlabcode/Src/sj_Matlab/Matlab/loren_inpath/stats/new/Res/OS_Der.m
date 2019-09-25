% OS_Der: Calculates derivatives for the David-Johnson approximation to the
%         variances and covariances of normal order statistics.  Based on AS 128.1,
%         Davis & Stephens (1978), translated from Fortran.  Called by OrderStats().
%

function [dx,d2x,d3x,d4x,d5x] = os_der(x)
  rad2pi = 2.506628274631d0;
  twopi = 2*pi;
     
  x2 = x*x;
  dx = rad2pi*exp(x2/2);
  d2x = twopi*x*exp(x2);
  d3x = twopi*rad2pi*(2*x2+1)*exp(1.5*x2);
  term = twopi*twopi*exp(2*x2);
  d4x = term*x*(6*x2+7);
  d5x = term*dx*(x2*(24*x2+46)+7);
      
  return;
      
