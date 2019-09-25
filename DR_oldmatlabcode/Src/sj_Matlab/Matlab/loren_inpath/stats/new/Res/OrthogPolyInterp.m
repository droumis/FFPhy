% OrthogPolyInterp: Interpret effect of varying orthogonal polynomial coefficients
%                   one at a time, via plots.
%
%     Usage: orthogpolyinterp(term,propdev,b,alpha,beta,x,y,{plotpts})
%
%         term =    identifier of coefficient to be varied (0-order).
%         propdev = proportional deviation below and above current value.
%         b =       column vector of regression coefficients.
%         alpha =   alpha coefficients used to calculate polynomials.
%         beta =    beta coefficients used to calculate polynomials.
%         x,y =     coordinates of original data.
%         plotpts = optional boolean flag indicating, if true, that
%                     original points are to be superimposed onto plots
%                     [default = 0].
%

% RE Strauss, 3/19/03

function orthogpolyinterp(term,propdev,b,alpha,beta,x,y,plotpts)
  if (~nargin) help orthogpolytinterp; return; end;
  
  if (nargin < 8) plotpts = []; end;
  
  if (isempty(plotpts)) plotpts = 0; end;
  
  term = term+1;
  xx = linspace(min(x),max(x))';
  bb = b;
  
  bb(term) = b(term)-propdev;
  ylow = orthogpolyregrp(xx,bb,alpha,beta);
  
  bb(term) = b(term)+propdev;
  yhigh = orthogpolyregrp(xx,bb,alpha,beta);
  
show  
  plot(xx,ylow,'k',xx,yhigh,'k');
  if (plotpts)
    hold on; 
    plot(x,y,'ko');
    hold off;
    x = x(:);
    y = y(:);
    putbnds([xx;xx;x],[ylow;yhigh;y]);
  else
    putbnds([xx;xx],[ylow;yhigh]);
  end;
  
  delta = 0.07;
  text(xx(end)+delta,ylow(end),'Low');
  text(xx(end)+delta,yhigh(end),'High');
  puttitle(sprintf('Coefficient %d',term-1));
  
  return;
  