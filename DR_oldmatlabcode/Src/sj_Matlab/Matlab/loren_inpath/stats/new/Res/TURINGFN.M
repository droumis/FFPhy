% TURINGFN: Rashevsky-Turing reaction-diffusion equations.  Defines the
%           differential equations to be integrated.  Returns the state
%           derivative vector, given the current time, state vector, and 
%           R-T constants.
%
%     Usage: yprime = turingfn(t,y,c)
%
%       t = time [not used in calculations of derivatives]
%       y = vector of state-variable values
%       c = Rashevsky-Turing constants
%

% RE Strauss, 6/4/94

function yprime = turingfn(t,y,c)

  K1 = c(1);      % Production of A
  K2 = c(2);      % Degradation of A
  K3 = c(3);      % Transformation from A to B
  K4 = c(4);      % Degradation of B
  K5 = c(5);      % Conc-independent production of A
                  %   (=Turing's "evocation parameter")
  K  = c(6);      % Saturation level of A
  D  = c(7);      % Diffusion coefficient for B

  ar = y(1);      % A, right side
  al = y(2);      % A, left side
  br = y(3);      % B, right side
  bl = y(4);      % B, left side

  dar = (K1-K3)*ar - K2*br*(ar/(ar+K)) + K5;
  dbr =      K3*ar - K4*br             + D*(bl-br);
  dal = (K1-K3)*al - K2*bl*(al/(al+K)) + K5;
  dbl =      K3*al - K4*bl             + D*(br-bl);

  yprime = zeros(4,1);             % Allocate return vector for
  yprime(1) = dar;                 %   state derivatives
  yprime(2) = dal;
  yprime(3) = dbr;
  yprime(4) = dbl;

  return;

