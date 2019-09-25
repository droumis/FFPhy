% TRAJFUN2: Provides the function to be minimized for the fitting of a
%     two-stage allometric trajectory function.
%
%     Usage:  [sse,b] = trajfun2(p,J,A)
%
%        p = parameter vector [s b1_j b1_a] where:
%                s = x,y-coordinates of inflection point for bilinear spline
%                b1_j = slope of juvenile portion of trajectory
%                b1_a = slope of adult portion of trajectory
%        J =   [Nj x 2] data matrix for juveniles
%        A =   [Na x 2] data matrix for adults
%        sse = total SSE as function of inflection point s
%        b =   [p x 4] vector of regression coefficients, in the form:
%                [b1_j, b0_j, b1_a, b0_a]
%
%     For J and A, the first column is a size measure and the second column
%     is the character to be fitted.
%

% RE Strauss, 6/13/93

function [sse,b] = trajfun2(p,J,A)
   X=1;                             % Coordinate indices
   Y=2;

   s = p(1:2);                      % Retrieve inflection point coords
   b1_j = p(3);                     % Retrieve juv slope
   b1_a = p(4);                     % Retrieve adult slope

   [rj,cj] = size(J);               % Matrix dimensions
   [ra,ca] = size(A);

   yj = J(:,2);                     % Dependent character var
   xj = [J(:,1) ones(rj,1)];        % Augmented independent size var
   ya = A(:,2);
   xa = [A(:,1) ones(ra,1)];

   bj = zeros(2,1);                 % Allocate coefficient vectors
   ba = zeros(2,1);

   bj(1) = b1_j;                    % Slopes
   ba(1) = b1_a;
   bj(2) = s(Y)-(bj(1).*s(X));      % Intercepts
   ba(2) = s(Y)-(ba(1).*s(X));
   b = [bj; ba]';                   % Vector to be returned

   ej = yj - xj*bj;                 % Residuals and SSE
   ea = ya - xa*ba;
   sse = ej'*ej + ea'*ea;

   % Assign a penalty to the SSE if out of the central box.
   % Penalty must be a function of the SSE to provide a gradient back to
   % the center.

   centr_j = mean(J);               % Group centroids
   centr_a = mean(A);
   d =  centr_a - centr_j;

   min_j = min(J);
   max_j = max(J);
   min_a = min(A);
   max_a = max(A);

   limit_left =  max_j(X);          % Modification for loricariid project
   limit_right = min_a(X);

%   limit_left =  centr_j(X) + 0.05*d(X);    % 5% inner horizontal buffer
%   limit_right = centr_a(X) - 0.05*d(X);
   limit_below = min_j(Y) - d(Y);           % 100% outer vertical buffer
   limit_above = max_a(Y) + d(Y);

   if (s(X)<limit_left  | s(X)>limit_right | ...
       s(Y)<limit_below | s(Y)>limit_above)
      sse = 100 * sse;             % Penalty
   end;

   return;

