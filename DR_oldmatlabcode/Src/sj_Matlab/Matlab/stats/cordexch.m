function [settings, X] = cordexch(nfactors,nruns,model,varargin)
%CORDEXCH D-Optimal design of experiments - coordinate exchange algorithm.
%   [SETTINGS, X] = CORDEXCH(NFACTORS,NRUNS,MODEL) generates a D-optimal
%   design having NRUNS runs for NFACTORS factors.  SETTINGS is the
%   matrix of factor settings for the design, and X is the matrix of
%   term values (often called the design matrix).  MODEL is an optional
%   argument that controls the order of the regression model.
%   MODEL can be any of the following strings:
%
%     'linear'        constant and linear terms (the default)
%     'interaction'   includes constant, linear, and cross product terms.
%     'quadratic'     interactions plus squared terms.
%     'purequadratic' includes constant, linear and squared terms.
%
%   Alternatively MODEL can be a matrix of term definitions as
%   accepted by the X2FX function.
%
%   [SETTINGS, X] = CORDEXCH(...,'PARAM1',VALUE1,'PARAM2',VALUE2,...)
%   provides more control over the design generation through a set of
%   parameter/value pairs.  Valid parameters are the following:
%
%      Parameter     Value
%      'bounds'      Lower and upper bounds for each factor, specified
%                    as a 2-by-NFACTORS matrix.  Alternatively, this value
%                    can be a cell array containing NFACTORS elements, each
%                    element specifying the vector of allowable values for
%                    the corresponding factor.
%      'categorical' Indices of categorical predictors.
%      'display'     Either 'on' or 'off' to control display of
%                    iteration counter (default = 'on').
%      'excludefun'  Function to exclude undesirable runs.
%      'init'        Initial design as an NRUNS-by-NFACTORS matrix
%                    (default is a randomly selected set of points).
%      'levels'      Vector of number of levels for each factor.
%      'maxiter'     Maximum number of iterations (default = 10).
%      'tries'       Number of times to try do generate a design from a
%                    new starting point, using random points for each
%                    try except possibly the first (default 1). 
%
%   The CORDEXCH function searches for a D-optimal design using a
%   coordinate exchange algorithm.  It creates a starting design, and then
%   iterates by changing each coordinate of each design point in an attempt
%   to reduce the variance of the coefficients that would be estimated
%   using this design.
%
%   If the 'excludefcn' function is F, it must support the syntax B=F(S) 
%   where S is a matrix of K-by-NFACTORS columns containing settings,
%   and B is a vector of K boolean values.  B(j) is true if the jth row
%   of S should be excluded.
%
%   Example:
%      % Design for two factors, quadratic model
%      sortrows(cordexch(2,9,'q'))
%
%      % Design for 2 of 3 factors making up a mixture, where factor
%      % values are up to 50%, and the two factors must not make up
%      % less than 15% or greater than 85% of the whole mixture
%      f = @(x) sum(x,2)>85 | sum(x,2)<15;
%      bnds = [0 0;50 50];
%      x=sortrows(cordexch(2,9,'q','bounds',bnds,'levels',101,'excl',f))
%      plot(x(:,1),x(:,2),'bo')
%
%   See also ROWEXCH, DAUGMENT, DCOVARY, X2FX.

%   The undocumented parameters 'start' and 'covariates' are used when
%   this function is called by the daugment and dcovary functions to
%   implement D-optimal designs with fixed rows or columns.
   
%   Copyright 1993-2005 The MathWorks, Inc. 
%   $Revision: 2.14.2.4 $  $Date: 2005/11/18 14:27:54 $

[eid,emsg,startdes,    covariates,   dodisp,     xinit,    maxiter,  ...
          tries,  bnds,         nlevels,    excluder, categ] = ...
                  doptargcheck('cordexch',{},nfactors,nruns, varargin{:});
if ~isempty(eid)
   error(sprintf('stats:cordexch:%s',eid),emsg);
end
quiet = isequal(dodisp,'off');

% Create design with starting rows, covariates, and zeros
nobs = size(startdes,1);
settings = zeros(nruns, nfactors);
if ~isempty(xinit)
   if ~isempty(excluder) && any(excluder(xinit))
       warning('stats:cordexch:ConstraintsViolated',...
               'Some points in the starting design violate the constraints.');
   end
end

% Add initial set of rows never to change, if any
if ~isempty(startdes)
   settings = [startdes; settings];
end

% Add fixed covariates, if any
if ~isempty(covariates)
   settings = [settings covariates];
end

if nargin == 2 || isempty(model)
   model = 'linear';
end
if ~ischar(model)
   if size(settings,2) ~= size(model,2)
      error('stats:cordexch:InputSizeMismatch',...
            'The number of columns in a numeric model matrix must equal the number of factors.');
   end
   modelorder = max(model(:));   % max order of any factor
else
   modelorder = 2;               % max of named models that we provide
end

if isempty(nlevels)
   nlevels = (modelorder+1) * ones(1,nfactors);
   nlevels(categ) = 2;
end
catlevels = nlevels(categ);
iscat = ismember(1:nfactors, categ);

% Convert from actual categorical factor levels to 1:maxlevel
if ~isempty(xinit)
    [xinit,eid,emsg] = levels2numbers(xinit,bnds,categ,nlevels);
    if ~isempty(eid)
        error(eid,emsg);
    end
end

% Create an array nxij containing all possible values for each factor
% if the number of levels is the same, or a cell array nxijcell if
% different factors have different numbers of levels
nxij = [];
nxijcell = {};
if all(nlevels==nlevels(1))
   usecell = false;
   nxij = zeros(nlevels(1),nfactors);
   for j=1:nfactors
      if iscell(bnds)
          bndsj = bnds{j};
      else
          bndsj = bnds(:,j);
      end
      if iscat(j)
          nxij(:,j) = (1:nlevels(j))';
      else
          nxij(:,j) = linspace(bndsj(1),bndsj(end),nlevels(j))';
      end
   end
   rowlist = zeros(nlevels(1),1);
else
   usecell = true;
   nxijcell = cell(nfactors,1);
   for j=1:nfactors
      if iscell(bnds)
          bndsj = bnds{j};
      else
          bndsj = bnds(:,j);
      end
      if iscat(j)
          nxijcell{j} = (1:nlevels(j))';
      else
          nxijcell{j} = linspace(bndsj(1),bndsj(end),nlevels(j))';
      end
   end
end

% Change exclusion function if necessary to deal with factor levels as numbers
if ~isempty(categ) && ~isempty(excluder)
    excluder = @(x) excluder(numbers2levels(x,bnds,categ,nlevels));
end

% Create Iteration Counter Figure.
if ~quiet
   [f,settry,setiter] = statdoptdisplay('Coordinate exchange');
end

% Generate designs, pick best one
bestdet = -Inf;
setstart = settings;
bestset = settings;
emsg = '';
eid = '';
warnedexcl = false;
warnedstart = false;
for j=1:tries
    % Update counter.
    if ~quiet
       settry(j);
    end

    % Create a new random starting design after the first iteration
    if isempty(xinit)
        [xinit,warnedexcl] = randstart(bnds, nlevels, nruns, nfactors, ...
                                       excluder, warnedexcl,  categ);
    end
    setstart(nobs+1:end,1:nfactors) = xinit;
    xinit = [];

    % Generate the term values
    [Xstart,model,termstart] = x2fx(setstart,model,categ,catlevels);

    % First time only, do error checking and gather some information
    if j==1
       totdf = size(Xstart,2);
       if totdf>(size(Xstart,1)+nobs)
          eid = 'stats:cordexch:TooFewRuns';
          emsg = 'There are not enough runs to fit the specified model.';
          break
       end
       modelorder = max(model,[],1);   % max order of each factor
       if any(~iscat & nlevels<modelorder(1:nfactors)+1)
          eid = 'stats:cordexch:TooFewLevels';
          emsg = 'Not all factors have enough levels to fit the model you specified.';
          break
       end

       temp = zeros(1,totdf);
       temp(termstart) = 1;
       termfactors = cumsum(temp);
       
       factorterms = false(nfactors,totdf);
       for row=1:nfactors
          factorterms(row,:) = ismember(termfactors,find(model(:,row)>0));
       end
    end

    % Generate a new design from here by coordinate exchange
    % (this nested function uses factorterms and some other variables)
    [newset,newX,logdetX] = gen1design(setstart,Xstart);
    if logdetX>bestdet
        bestdet = logdetX;
        bestset = newset;
        X = newX;
    end
end

if ~quiet && ishandle(f)
   close(f);
end

if ~isempty(eid)
   error(eid,emsg);
end

% Convert to actual categorical factor levels
settings = numbers2levels(bestset,bnds,categ,nlevels);


% -------------------
function [settings,X,logdetX] = gen1design(settings,X)

% This is a nested function that shares some variables with its caller

[Q,R]=qr(X,0);

% Adjust starting design if it is rank deficient, because the algorithm
% will not proceed otherwise.
if rank(R)<size(R,2)
   if ~warnedstart
      warning('stats:cordexch:BadStartingDesign',...
              'Starting design is rank deficient');
      warnedstart = true;
   end
   R = adjustr(R);
   wasbad = 1;
else
   wasbad = 0;
end
logdetX = 2*sum(log(abs(diag(R))));

iter = 0;
madeswitch = 1;
dcutoff = 1 + sqrt(eps);

while madeswitch > 0 && iter < maxiter
   madeswitch = 0;
   iter = iter + 1;
   
   % Update iteration counter.
   if ~quiet
      setiter(iter);
   end
   
   %Loop over rows of factor settings matrix.
   collist = randperm(nfactors);
   for row = (nobs+1):(nobs+nruns)
      fx = X(row,:);
      E = [];
      %Loop over columns of factor settings matrix.
      for col = collist
         if usecell
            newset = nxijcell{col};
            xnew = repmat(settings(row,:),numel(newset),1);
            xnew(:,col) = newset;
         else
            rowlist(:) = row;
            xnew = settings(rowlist,:);
            xnew(:,col) = nxij(:,col);
         end

         if ~isempty(excluder)
            excluderows = excluder(xnew);
            xnew = xnew(~excluderows,:);
            if isempty(xnew)
                continue
            end
         end

         % Update affected terms only
         t = find(model(:,col)>0);
         fxnew = fx(ones(size(xnew,1),1),:);
         fxnew(:,factorterms(col,:)) = x2fx(xnew,model(t,:),categ,catlevels);

         % Compute change in determinant.
         if isempty(E)
            E = fx/R;
            dxold = E*E';
         end
         F = fxnew/R;
         dxnew = sum(F.*F,2);
         dxno  = F*E';

         d = (1 + dxnew).*(1 - dxold) + dxno.^2;

         % Find the maximum change in the determinant, switch if >1
         [d,idx] = max(d);
         if d>dcutoff || (iter==1 && ~iscat(col))
            madeswitch = 1;
            logdetX = log(d) + logdetX;
            settings(row,col) = xnew(idx,col);
            X(row,:) = fxnew(idx,:);
            fx = X(row,:);
            [Q,R] = qr(X,0);
            if wasbad
               if rank(R)<size(R,2)
                  R = adjustr(R);
               else
                  wasbad = 0;
               end
            end
            E = [];
         end
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

% --------------------------------------------------
function [xinit,warnedexcl] = randstart(bnds, nlevels, nruns, nfactors, ...
                                        excluder, warnedexcl, categ)
%RANDSTART Create random but feasible initial design
xinit = zeros(0,nfactors);
startiter = 1;
if isempty(excluder)
   blocksize = nruns;
else
   blocksize = max(1000,nruns);
end
catcols = ismember(1:nfactors, categ);
contcols = ~catcols;
ncatcols = sum(catcols);
ncontfactors = sum(contcols);
if iscell(bnds)  % only need extremes here, convert to matrix for convenience
   tempbnds = zeros(2,nfactors);
   for j=1:nfactors
       v = bnds{j};
       tempbnds(1,j) = v(1);
       tempbnds(2,j) = v(end);
   end
   bnds = tempbnds;
end
bnds(1,catcols) = 1;
bnds(2,catcols) = nlevels(catcols);

% Repeatedly attempt to generate points, discarding excluded ones
while(size(xinit,1)<nruns && startiter<=20)
    block = zeros(blocksize,nfactors);
    block(:,contcols) = repmat(bnds(1,contcols),blocksize,1) + ...
            repmat(diff(bnds(:,contcols),1,1),blocksize,1) ...
                .* rand(blocksize,ncontfactors);
    block(:,catcols) = ceil(rand(blocksize,ncatcols) ...
                            .* repmat(bnds(2,catcols),blocksize,1));
    if ~isempty(excluder)
       excluderows = excluder(block);
       if startiter==20 % last time only, save discards in case needed
           discards = block(excluderows,:);
       end
       block(excluderows,:) = [];
    end       
    xinit = [xinit; block];
    startiter = startiter + 1;
end
if size(xinit,1)<nruns
    % Use bad points if absolutely necessary, hoping they can be removed
    % during the coordinate exchange process
    if ~warnedexcl
       warning('stats:cordexch:ConstraintsViolated',...
             'Cannot find random initial design to satisfy the constraints.');
       warnedexcl = true;
    end
    N = nruns - size(xinit,1);
    xinit = [xinit; discards(1:N,:)];
else
    xinit = xinit(1:nruns,:);
end
end

% -------------------------------------------------------
function [xinit,eid,emsg] = levels2numbers(xinit,bnds,categ,nlevels)
% Renumber levels from values in bnds to numbers in 1:nlevels

eid = '';
emsg = '';
    
for j=1:length(categ)
    cj = categ(j);
    if iscell(bnds)
       setlist = bnds{cj};
    else
       setlist = linspace(bnds(1,cj),bnds(2,cj),nlevels(cj));
    end
    [foundit,catrow] = ismember(xinit(:,cj),setlist);
    if any(~foundit)
       eid = 'stats:cordexch:BadInit';
       emsg = sprintf(...
           ['One or more entries in column %d of the initial design do not'...
            '\nappear in the list specified by the ''bounds'' parameter.'],...
             cj);
       return
    end
    xinit(:,cj) = catrow;
end
end

% -------------------------------------------------------
function settings = numbers2levels(settings,bnds,categ,nlevels)
% Renumber levels from numbers in 1:nlevels to values in bnds

for j=1:length(categ)
    cj = categ(j);
    if iscell(bnds)
        setlist = bnds{cj};
    else
        setlist = linspace(bnds(1,cj),bnds(2,cj),nlevels(cj));
    end
    settings(:,cj) = setlist(settings(:,cj));
end
end