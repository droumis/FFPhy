% ODE:  Modification of the MATLAB function ODE45 to solve differential 
%       equations using the higher-order method (4th and 5th order 
%       Runge-Kutta formulas).  Modified to permit the passing of a vector 
%       of arguments to the user-function defining the differential equations.
%
%     Usage: [tout,yout] = ode(ypfunc, t0, tfinal, y0, c)
%
%       ypfunc - string containing name of user-supplied problem description
%                 Call: yprime = ypfunc(t,y)
%                    ypfunc - 'funcname'
%                    t      - time (scalar)
%                    y      - solution column-vector
%                    yprime - returned derivative column-vector:
%                               yprime(i) = dy(i)/dt.
%       t0     - initial value of t
%       tfinal - final value of t
%       y0     - initial value column-vector (one per equation)
%       c      - vector of constants passed to ypfunc for solution of 
%                  differential equations
%       tout   - returned integration time points (column-vector)
%       yout   - returned solution, one solution column-vector per tout-value
%
%       The result can be displayed by: plot(tout, yout).
%

% RE Strauss, 4/20/95

function [tout,yout] = ode(ypfun,t0,tfinal,y0,c)
  tol = 1.e-6;
  trace = 0;                          % Set to 1 to print intermediate results

  % The Fehlberg coefficients:
  alpha = [1/4  3/8  12/13  1  1/2]';
  beta  = [ [    1      0      0     0      0    0]/4
            [    3      9      0     0      0    0]/32
            [ 1932  -7200   7296     0      0    0]/2197
            [ 8341 -32832  29440  -845      0    0]/4104
            [-6080  41040 -28352  9295  -5643    0]/20520 ]';
  gamma = [ [902880  0  3953664  3855735  -1371249  277020]/7618050
            [ -2090  0    22528    21970    -15048  -27360]/752400 ]';
  pow = 1/5;

  % Initialization
  t = t0;
  hmax = (tfinal - t)/16;
  h = hmax/8;
  y = y0(:);
  f = zeros(length(y),6);
  chunk = 128;
  tout = zeros(chunk,1);
  yout = zeros(chunk,length(y));
  k = 1;
  tout(k) = t;
  yout(k,:) = y.';

  if (trace)
    clc, t, h, y
  end;

  % The main loop
  while (t < tfinal) & (t + h > t)
    if (t + h > tfinal)
      h = tfinal - t;
    end

    % Compute the slopes
    temp = feval(ypfun,t,y,c);
    f(:,1) = temp(:);
    for j = 1:5
      temp = feval(ypfun, t+alpha(j)*h, y+h*f*beta(:,j), c);
      f(:,j+1) = temp(:);
    end;

    % Estimate the error and the acceptable error
    delta = norm(h*f*gamma(:,2),'inf');
    tau = tol*max(norm(y,'inf'),1.0);

    % Update the solution only if the error is acceptable
    if delta <= tau
      t = t + h;
      y = y + h*f*gamma(:,1);
      k = k+1;
      if k > length(tout)
        tout = [tout; zeros(chunk,1)];
        yout = [yout; zeros(chunk,length(y))];
      end
      tout(k) = t;
      yout(k,:) = y.';
    end;
    if (trace)
      home, t, h, y
    end;

    % Update the step size
    if delta ~= 0.0
      h = min(hmax, 0.8*h*(tau/delta)^pow);
    end;
  end;

  if (t < tfinal)
    disp('Singularity likely.')
    t
  end

  tout = tout(1:k);
  yout = yout(1:k,:);
  return;
