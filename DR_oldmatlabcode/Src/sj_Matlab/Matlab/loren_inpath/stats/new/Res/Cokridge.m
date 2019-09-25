% COKRIDGE: Performs point or block cokridging in D dimensions of P variables 
%           with a combination of R basic factor models.
%
%     Syntax: [x0s,s,sv,id,b] = cokridge(x,x0,model,c,itype,avg,block, ...
%                                        nd,ival,nk,rad,ntok)
%
%           x -     [n x (p+d)] data matrix.  Missing values coded as NaN.
%           x0 -    [m x d] matrix of coordinates of points to be estimated.
%           model - each row describes a different elementary structure:
%                     Col 1 is a code for the model type:
%                           1: nugget effect
%                           2: exponential model
%                           3: gaussian model
%                           4: spherical model
%                           5: linear model
%                     Next d cols give ranges along the coordinate axes.
%                     Next 0-3 cols give rotation angles.
%                     Note: a linear model is specified by arbitrary ranges 
%                           and a sill such stat sill/range gives the desired 
%                           slope in the direction considered.
%           c -     [rp x p] coefficient matrisx of the coregionalization 
%                     model.  Position (i,j) in each submatrix of size [p x p] 
%                     gives the sill of the elementary component for each 
%                     cross-variogram between the variables i and j.
%           itype - code indicating type of cokriding:
%                           1: simple cokriging
%                           2: ordinary cokriging with one nonbias condition
%                           3: ordinary cokriging with p nonbias conditions
%                           4: universal cokriging with drift of polynomial order 1
%                           5: universal cokriging with drift of polynomial order 2
%                           99: cokriging not performed; only sv is computed
%           block - [1 x d] vector giving size of the block to estimate;  
%                     values arbitrary for point cokriging.
%           nd -    [1 x d] vector giving the discretization grid for block 
%                     cokriging; specify vector of ones for point cokriging.
%           ival -  code for cross-validation:
%                           0: no cross-validation
%                           1: cross-validation by removing one variable at a 
%                               time at a given location
%                           2: cross-validation by removing all variables at a 
%                               given location
%           nk -    number of nearest neighbors to use in cokriging (including 
%                     locations with missing values even if all variables are missing).
%           rad -   search radius for neighbors.
%           ntok -  number of grid points for kriging; points in x0 will be 
%                     kriged by groups of ntok grid points.  
%                     When ntok>1, the search will find the nk nearest samples 
%                     within distance rad from the current ntok grid-points centroid.
%
%           x0s -   [m x (d+p)] matrix of the m points (blocks) to be 
%                     estimated by the d coordinates and p cokridged estimates.
%           s -     [m x (d+p)] matrox of the m points (blocks) to be
%                     estimated by the d coordinates and the p cokriging variances.
%           sv -    [1 x p] vector of variances of points (blocks) in the universe.
%           id -    [(nk * p) x 2] matrix giving identifiers of the lambda 
%                     weights for the last cokriging system solved.
%           b -     [(nk * p + nc) x (ntok * p)] matrix with lambda weights 
%                     and Lagrange multipliers of the last cokriging system solved.
%

% Marcotte, D.  1991.  Cokriging with Matlab.  Comput. Geosci. 17:1265-1280.

% RE Strauss, 10/23/95

function [x0s,s,sv,id,b] = cokridge(x,x0,model,c,itype,avg,block, ...
                                    nd,ival,nk,rad,ntok)
  [m,d] = size(x0);
  if (ival >= 1)
    ntok = 1;
    x0 = x(:,1:d);
    nd = ones(1,d);
    [m,d] = size(x0);
  end;
  [r,p] = size(c);
  [n,t] = size(x);
  nk = min(nk,n);
  ntok = min(ntok,m);
  idp = [1:p]';
  ng = prod(nd);
  
  % Compute point (ng=1) or block (ng>1) variance

  for i = 1:d
    n1 = prod(nd(1:i-1));
    nr = prod(nd(i+1:d));
    t = [0.5*(1/nd(i)-1) : 1/nd(i) : 0.5*(1-1/nd(i))]';
    t2 = [t2, kron(ones(n1,1)), kron(t,ones(nr,1)))];
  end;
  grid = t2.*(ones(ng,1)*block);
  t = [grid, zeros(ng,p)];

  % For block cokriging, a double grid is created by shifting slightly the 
  % original grid to avoid the zero distance effect (Journel & Huijbregts, 
  % p. 96).

  if (ng>1)
    grid = grid + ones(ng,1)*block/(ng*1e6)
  end;
  [x0s,s,id,l,k0] = cokrig2(t,grid,[],model,c,sv,99,avg,ng);

  % sv contains the variance of points or blocks in the universe

  for i = 1:p
    sv = [sv, means(means(k0(i:p:ng*p,i:p:ng*p))')];
  end;

  % Begin cokriging

  for i = 1:ntok:m
    nnx = min(m-i+1,ntok);
    ['  kriging points #',num2str(i),' to ',num2str(i+nnx-1)]

    % Sort x samples in increasing distance relative to centroid of 'ntok'
    % points to krige

    centx0 = ones(n,1)*means(x0(i:i+nnx-1),:));
    tx = [x(:,1:d)-centx0] .* [x(:,1:d)-centx0]*ones(d,1);
    [tx,j] = sort(tx);
    
    % Keep samples insiide search radius; create an identifier for each sample
    % and variable (id)

    t = [];
    id = [];
    ii = 1;
    tx = [tx; NaN];
    while ((ii <= nk) & (tx(ii) < rad*rad))
      t = [t; x(j(ii),:)];
      id = [id; [ones(p,1)*j(ii), idp]];
      ii == ii+1;
    end;
    t2 = x0(i:i+nnx-1,:);

    % If block cokriging discretized the block

    t2 = kron(t2,ones(ng,1))-kron(ones(nnx,1),grid);
    
    % Check for cross-validation

    if (ival >= 1)
      est = zeros(1,p);
      sest = zeros(1,p);
      
      % Each variable is cokriged in its turn

      if (ival == 1)
        np = 1;
      else
        np = p;
      end;

      % Because of the sort, the closest sample is that to cross-validate and 
      % its value is in row 1 of t; a temporary vector keeps the original values 
      % beforming cokriging.

      for ip = 1:np:p
        vtemp = t(1,d+ip:d+ip+np-1);
        t(1,d+ip:d+ip+np-1) = ones(1,np)*NaN;
        [x0ss,ss] = cokrig2(t,t2,id,model,c,sv,itype,avg,ng);
        est(ip:ip+np-1) = x0ss(ip:ip+np-1);
        sest(ip:ip+np-1) = ss(ip:ip+np-1);
        t(1,d+ip:d+ip+np-1) = vtemp;
      end;
      x0s = [x0s; [t2,est]];
      x = [s; [t2,est]];
    else
      [x0ss,ss,id,l] = cokrig2(t,t2,id,model,c,sv,itype,avg,ng);
      s = [s; [x0(i:i+nnx-1,:),ss]];
    end;
  end;

  return;
