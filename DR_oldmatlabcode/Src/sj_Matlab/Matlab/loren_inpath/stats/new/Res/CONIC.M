% CONIC: Plots the 6-parameter form of an arbitrary conic:
%           Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
%
%        Conservative upper and lower bounds on x can be specified,
%        and will be adjusted inward.
%
%     Usage: conic(A,B,C,D,E,F,{lower},{upper})
%
%             A,B,C,D,E,F = conic coefficients
%             lower =       lower bound for x [default =  50]
%             upper =       upper bound for x [default = -50]
%

% RE Strauss, 9/20/97
%   9/7/99 -  changed plot colors for Matlab v5.
%   1/4/00 -  changed usage of sqplot().

function conic(A,B,C,D,E,F,lower,upper)
  if (nargin < 7)
    lower = -50;
  end;
  if (nargin < 8)
    upper = 50;
  end;

  x = linspace(lower,upper,5000)';

  a = C*ones(size(x));          % Quadratic coefficients
  b = B*x + E;
  c = A*x.*x + D*x + F;

  r = 0.5./a;                   % Solve quadratic
  s = sqrt(b.*b-4.*a.*c);
  ypos = r.*(-b+s);
  yneg = r.*(-b-s);

  i = find(~imag(ypos) | ~imag(yneg));    % Find bounds of soln
  leni = length(i);

  if (leni == 0)
    disp('  No solution found');
  else
    lower = x(max([min(i),1]));
    upper = x(min([max(i),length(x)]));
    x = linspace(lower,upper,200)';  % New x bounds

    a = C*ones(size(x));          % Quadratic coefficients for y
    b = B*x + E;
    c = A*x.*x + D*x + F;

    r = 0.5./a;                   % Solve quadratic for y as fn(x)
    s = sqrt(b.*b-4.*a.*c);
    ypos = r.*(-b+s);
    yneg = r.*(-b-s);
    iypos = ~imag(ypos);
    iyneg = ~imag(yneg);

    upper = max(ypos(iypos));
    lower = min(yneg(iyneg));
    y = linspace(lower,upper,200)'; % y bounds

    a = A*ones(size(y));          % Quadratic coefficients for x
    b = B*y + D;
    c = C*y.*y + E*y + F;

    r = 0.5./a;                   % Solve quadratic for x as fn(y)
    s = sqrt(b.*b-4.*a.*c);
    xpos = r.*(-b+s);
    xneg = r.*(-b-s);
    ixpos = ~imag(xpos);
    ixneg = ~imag(xneg);

    close all;
    plot(x(iypos),ypos(iypos),'k');       % Plot pos y as fn(x)
    config = [x(iypos) ypos(iypos)];
    hold on;

    plot(x(iyneg),yneg(iyneg),'k');       % Plot neg y as fn(x)
    config = [config; x(iyneg) yneg(iyneg)];

    plot(xpos(ixpos),y(ixpos),'k');       % Plot pos x as fn(y)
    config = [config; xpos(ixpos) y(ixpos)];

    plot(xneg(ixneg),y(ixneg),'k');       % Plot neg x as fn(y)
    config = [config; xneg(ixneg) y(ixneg)];

    sqplot(config);                       % Square plot, equal axes
    hold off;
  end;

  return;
