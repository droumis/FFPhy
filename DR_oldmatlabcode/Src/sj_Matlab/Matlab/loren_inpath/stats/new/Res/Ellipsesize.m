% EllipseSize: area and perimeter of an ellipse
%
%     Usage:  [area,perim] = ellipsesize(a,b)
%
%         a,b =   vectors of min and max radii of ellipse (half-lengths of 
%                   major and minor axes, respectively).
%         ----------------------------------------------------------------
%         area =  area of ellipse.
%         perim = perimeter of ellipse.
%

% Ref: Beta Mathematics Handbook (2nd ed), pp. 80, 252.

% RE Strauss, 12/12/01

function [area,perim] = ellipsesize(a,b)
  get_perim = 0;
  if (nargout>1)
    get_perim = 1;
  end;

  ax = a(:);
  bx = b(:);
  if (length(a)~=length(b))
    error('  ELLIPSESIZE: input vectors must be same length');
  end;
  
  a = max([ax';bx'])';
  b = min([ax';bx'])';

  area = a.*b.*pi;
  
  if (get_perim)
    perim = zeros(size(area));
    k = sqrt((a.*a)-(b.*b))./a;
    for i = 1:length(k)                         % Elliptic integral
      perim(i) = 4*a(i)*quad('Ellipsesizef',0,pi/2,[],[],k);
    end;
  end;

  return;