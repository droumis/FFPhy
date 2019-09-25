% MNBRAK: Brackets a minimum for one-dimensional optimization.  Given an 
%         objective-function 'func', and distinct initial points ax and bx, 
%         this routine searches in the downhill direction (defined by the 
%         function as evaluated at the initial points) and returns new points 
%         [ax bx cx] that brack a minimum of the function.  Also returned are 
%         the function values at the three points, [fa fb fc].
%

% RE Strauss, 4/19/95

function [x,fx] = mnbrak(evalstr,x)
  gold = 1.618034;                    % Default ratio by which successive
                                      %   intervals are magnified
  glimit = 100;                       % Max mag for parabolic-fit step

  ax = x(1);
  bx = x(2);

  p = ax;
  fa = eval(evalstr);
  p = bx;
  fb = eval(evalstr);

  if (fb > fa)                        % Will go downhill from ax to bx,
    dum = ax;                         %   so switch points if necessary
    ax = bx;
    bx = dum;
    dum = fa;
    fa = fb;
    fb = dum;
  end;
  cx = bx + gold*(bx-ax);             % First guess at cx
  p = cx;
  fc = eval(evalstr);

  while (fb > fc)                     % Continue until bracket
    r = (bx-ax)*(fb-fc);              % Compute u by parabolic extrapolation
    q = (bx-cx)*(fb-fa);
    u = bx-(((bx-cx)*q-(bx-ax)*r)/(2*sign(q-r)*max(abs(q-r),eps)));

  end; % while (fb > fc)

  return;
