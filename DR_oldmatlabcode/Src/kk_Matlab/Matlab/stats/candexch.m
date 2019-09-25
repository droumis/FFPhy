function rowlist = candexch(fxcand,nruns,varargin)
%CANDEXCH D-optimal design from candidate set using row exchanges.
%   RLIST = CANDEXCH(C,NROWS) uses a row-exchange algorithm to select a
%   D-optimal design from the candidate set C.  C is an N-by-P matrix
%   containing the values of P model terms at each of N points.  NROWS
%   is the desired number of rows in the design.  RLIST is a vector of
%   length NROWS listing the selected rows.
%
%   The CANDEXCH function selects a starting design X at random, and
%   uses a row-exchange algorithm to iteratively replace rows of X by
%   rows of C in an attempt to improve the determinant of X'*X.
%
%   RLIST = CANDEXCH(C,NROWS,'PARAM1',VALUE1,'PARAM2',VALUE2,...)
%   provides more control over the design generation through a set of
%   parameter/value pairs.  Valid parameters are the following:
%
%      Parameter    Value
%      'display'    Either 'on' or 'off' to control display of
%                   iteration number (default = 'on').
%      'init'       Initial design as an NROWS-by-P matrix (default
%                   is a random subset of the rows of C).
%      'maxiter'    Maximum number of iterations (default = 10).
%      'tries'      Number of times to try do generate a design from a
%                   new starting point, using random points for each
%                   try except possibly the first (default 1). 
%
%   The ROWEXCH function also generates D-optimal designs using a
%   row-exchange algorithm, but it automatically generates a candidate
%   set that is appropriate for a specified model.
%
%   Example:  generate a D-optimal design when there is a restriction
%   on the candidate set, so the ROWEXCH function isn't appropriate.
%      F = (fullfact([5 5 5])-1)/4;   % factor settings in unit cube
%      T = sum(F,2)<=1.51;            % find rows matching a restriction
%      F = F(T,:);                    % take only those rows
%      C = [ones(size(F,1),1) F F.^2];% compute model terms including
%                                     % a constant and all squared terms
%      R = candexch(C,12);            % find a D-optimal 12-point subset
%      X = F(R,:);                    % get factor settings
%
%   See also CANDGEN, ROWEXCH, CORDEXCH, X2FX.

%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 1.3.4.4 $  $Date: 2005/11/18 14:27:47 $

maxiter = 10;
iter = 0;
madeswitch = 1;
dcutoff = 1 + sqrt(eps(class(fxcand)));

pnames = {'display'   'init'    'maxiter'    'tries'};
[eid,emsg,dodisp,xinit,maxiter, tries] = ...
             doptargcheck('candexch',pnames,size(fxcand,2),nruns,varargin{:});
if ~isempty(eid)
   error(eid,emsg);
end

% Create Iteration Counter Figure.
quiet = isequal(dodisp,'off');
if ~quiet
   [f,settry,setiter] = statdoptdisplay('Row exchange');
end

% Generate designs, pick best one
bestdet = -Inf;
rowlist = [];
for j=1:tries
    % Update counter.
    if ~quiet
       settry(j);
    end

    if isempty(xinit)
       X = fxcand(unidrnd(size(fxcand,1),nruns,1),:);
    else
       X = xinit;
       xinit = [];
    end
    [rows,logdetX] = gen1design(X);
    if logdetX>bestdet
        bestdet = logdetX;
        rowlist = rows;
    end
end

if ~quiet && ishandle(f)
    close(f);
end

% -------------------
function [rowlist,logdetX] = gen1design(X)

rowlist = zeros(nruns,1);
iter = 0;
madeswitch = 1;

F = [];
[Q,R]=qr(X,0);

% Adjust starting design if it is rank deficient, because the algorithm
% will not proceed otherwise.
if rank(R)<size(R,2)
   warning('stats:cordexch:BadStartingDesign',...
           'Starting design is rank deficient');
   R = adjustr(R);
   wasbad = 1;
else
   wasbad = 0;
end
logdetX = 2*sum(log(abs(diag(R))));

while madeswitch > 0 && iter < maxiter
   madeswitch = 0;
   iter = iter + 1;

   % Update iteration counter.
   if ~quiet
      setiter(iter);
   end
  
   for row = 1:nruns
      % Compute determinant factor over whole candidate set if not done yet
      if isempty(F)
         F = fxcand/R;
         dxcand = sum(F.*F, 2);
      end
      
      E = X(row,:)/R;
      dxold = E*E';
      dxswap  = F*E';
      
      dd = (1+dxcand) * (1-dxold) + dxswap.^2;
     
      % Find the maximum change in the determinant.
      [d,idx] = max(dd);
     
      % Switch rows if the maximum change is greater than 1.
      if (d > dcutoff) || (rowlist(row) == 0)
         madeswitch = 1;
         logdetX = log(d) + logdetX;
         rowlist(row) = idx;
         X(row,:) = fxcand(idx,:);
         [Q,R] = qr(X,0);
         if wasbad
            if rank(R)<size(R,2)
               R = adjustr(R);
            else
               wasbad = 0;
            end
         end
         F = [];  % needs re-computing using new R
      end
   end
end

end % of nested function
end % of outer function

% --------------------------------------------
function R = adjustr(R)
%ADJUSTR Adjust R a little bit so it will be non-singular

diagr = abs(diag(R));
p = size(R,2);
smallval = sqrt(eps)*max(diagr);
t = (diagr < smallval);
if any(t)
   tind = (1:p+1:p^2);
   R(tind(t)) = smallval;
end
end

