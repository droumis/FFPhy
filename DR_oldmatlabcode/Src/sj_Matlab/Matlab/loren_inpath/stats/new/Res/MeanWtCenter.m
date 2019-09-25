% MeanWtCenter: Calculates the weighted mean of a variable (optionally by group) where the weights
%               are the reciprocals of absolute deviations from the unweighted mean of a variable y,
%               raised to a power 'a'. For the default of a==1, the weights are the reciprocal 
%               distances from the mean of y.
%                 Example of use: calculate weighted means of a set of variables based on distances
%               from mean principal-component or discriminant-function scores.
%
%     Usage: [m,grpid] = meanwtcenter(X,y,{g},{a})
%
%         X =     [n x p] data matrix.
%         y =     corresonding [n x 1] vector of scores.
%         g =     optional [n x 1] group-membership vector for k groups.
%         a =     optional power to which absolute deviations from the mean are to be 
%                   raised to calculate weights [default = 1].
%         ---------------------------------------------------------------------------
%         m =     [k x p] matrix of weighted means.
%         grpid = [k x 1] vector of corresponding group identifiers.
%

% RE Strauss, 3/21/03

function [m,grpid] = meanwtcenter(X,y,g,a)
  if (~nargin) help meanwtcenter; return; end;
  
  if (nargin < 3) g = []; end;
  if (nargin < 4) a = []; end;
  
  [N,P] = size(X);
  
  if (isempty(g)) g = ones(N,1); end;
  if (isempty(a)) a = 1; end;
    
  ug = uniquef(g);
  ngrps = length(ug);
  grpid = ug;
  
  m = zeros(ngrps,P);
  
  for ig = 1:ngrps
    for ip = 1:P
      xi = X(g==ug(ig),ip);
      yi = y(g==ug(ig));
      w = 1./(abs(yi-mean(yi)).^a);
      m(ig,ip) = meanwt(xi,w);
    end;
  end;
  
  return;
  