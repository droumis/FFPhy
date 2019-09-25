function [Param,Covar,Resid,VarParam,Objective] = mvregress(Design, Data, varargin)
%MVREGRESS Multivariate regression with missing data.
%   [BETA,SIGMA,RESID]=MVREGRESS(X,Y) performs multivariate regression of
%   the multivariate observations in the N-by-D matrix Y on the predictor
%   variables in X, and returns a P-by-1 column vector BETA of coefficient
%   estimates, a D-by-D matrix SIGMA of the estimated covariance of Y, and an
%   N-by-D matrix RESID of residuals.  NaN values in X or Y are taken to be
%   missing.  Observations with missing values in X are ignored.  Missing
%   values in Y are handled according to the value of the 'algorithm'
%   parameter described below.
%
%   Y is an N-by-D matrix of D-dimensional multivariate observations.  X may
%   be either a matrix or a cell array.  If D=1, X may be an N-by-P design
%   matrix of predictor variables.  For any value of D, X may also be a cell
%   array of length N, each cell containing a D-by-P design matrix for one
%   multivariate observation.  If all observations have the same D-by-P
%   design matrix, X may be a single cell.
%
%   [BETA,SIGMA,RESID,VARPARAM]=MVREGRESS(...) also returns an estimated
%   covariance matrix of the estimates.  By default, or if the 'varformat'
%   parameter is 'beta' (see below), VARPARM is the estimated covariance
%   matrix of the coefficient estimates BETA.  If the 'varformat' parameter
%   is 'full', VARPARAM is the estimated covariance matrix of the combined
%   BETA and SIGMA estimates.
%
%   [BETA,SIGMA,RESID,VARPARAM,OBJECTIVE]=MVREGRESS(...) also returns the
%   value of the objective function, or log likelihood, after the last
%   iteration.
%
%   [...]= MVREGRESS(X,Y,'PARAM1',VALUE1,'PARAM2',VALUE2,...) specifies
%   additional parameter name/value pairs chosen from the following:
%     'algorithm'  Either 'ecm' to compute the maximum likelihood estimates
%                  via the ECM algorithm, 'cwls' to perform least squares
%                  (optionally conditionally weighted by an input covariance
%                  matrix), or 'mvn' to omit observations with missing data
%                  and compute the ordinary multivariate normal estimates.
%                  Default is 'mvn' for complete data, 'ecm' for missing data
%                  when the sample size is sufficient to estimate all
%                  parameters, and 'cwls' otherwise.
%     'maxiter'    Maximum number of iterations (default 100).
%     'tolbeta'    Convergence tolerance for BETA (default sqrt(eps)).
%                  Iterations continue until the TOLBETA and TOLOBJ conditions
%                  are met.  The test for convergence at iteration k is
%                     ||BETA(k)-BETA(k-1)|| < sqrt(P)*TOLBETA * (1+||BETA(k)||)
%                  where ||v|| represents the norm of the vector v.
%     'tolobj'     Convergence tolerance for changes in the objective function
%                  (default eps^(3/4)).  The test is
%                     |Obj(k)-Obj(k-1)| < TolObj * (1 + |Obj(k)|)
%                  If both TOLOBJ and TOLBETA are 0, the function performs
%                  MAXITER iterations with no convergence test.
%     'param0'     A vector of P elements to be used as the initial estimate
%                  for PARAM.  Default is a zero vector.  Not used for the
%                  'mvn' algorithm.
%     'covar0'     A D-by-D matrix to be used as the initial estimate for
%                  SIGMA.  Default is the identity matrix.  For the 'cwls'
%                  algorithm, this matrix is usually a diagonal matrix and it
%                  is not changed during the iterations, so the input value
%                  is used as the weighting matrix at each iteration.
%     'outputfcn'  An output function.
%     'varformat'  Either 'beta' to compute VARPARAM for BETA only (default),
%                  or 'full' to compute VARPARAM for both BETA and SIGMA.
%     'vartype'    Either 'hessian' to compute VARPARAM using the Hessian or
%                  observed information (default), or 'fisher' to compute the
%                  complete-data Fisher or expected information.  The 'hessian'
%                  method takes into account the increased uncertainties due
%                  to missing data, while the 'fisher' method does not.
%
%   The RESID values corresponding to missing values in Y are the differences
%   between the conditionally-imputed values for Y and the fitted values.  The
%   SIGMA estimate is not the sample covariance matrix of the RESID matrix.
%
%   The output function is called with three arguments:
%      1.  Vector of current parameter estimates
%      2.  A structure with fields 'Covar' for the current value of the
%          covariance matrix, 'iteration' for the current iteration number,
%          and 'fval' for the current value of the objective function.
%      3.  A text string that is 'init' when called during initialization,
%          'iter' when called after an iteration, and 'done' when called
%          after completion.
%
%   See also REGSTATS, MANOVA1, MVREGRESSLIKE.

% References:
%    [1] Roderick J. A. Little and Donald B. Rubin, Statistical Analysis with
%        Missing Data, 2nd ed., John Wiley & Sons, Inc., 2002.
%    [2] Xiao-Li Meng and Donald B. Rubin, "Maximum Likelihood Estimation via
%        the ECM Algorithm," Biometrika, Vol. 80, No. 2, 1993, pp. 267-278.
%    [3] Joe Sexton and Anders Rygh Swensen, "ECM Algorithms that Converge at
%        the Rate of EM," Biometrika, Vol. 87, No. 3, 2000, pp. 651-662.
%    [4] A. P. Dempster, N.M. Laird, and D. B. Rubin, "Maximum Likelihood from
%        Incomplete Data via the EM Algorithm," Journal of the Royal Statistical
%        Society, Series B, Vol. 39, No. 1, 1977, pp. 1-37.

%    Copyright 2006-2007 The MathWorks, Inc.
%    $Revision: 1.1.6.4.2.1 $ $Date: 2007/01/30 02:22:24 $

% Step 1 - check arguments
if nargin < 2
    error('stats:mvregress:MissingInputArg', ...
          'Missing required arguments X or Y.');
end
if isempty(Data)
    error('stats:statcheckmvnr:EmptyDataArray', ...
          'Y array is empty - cannot continue.');
end

if isempty(Design)
    error('stats:statcheckmvnr:EmptyDesignArray', ...
          'X array is empty - cannot continue.');
end

okargs =   {'maxiter'   'tolparam'  'tolobj'    'param0'  'covar0' ...
            'algorithm' 'outputfcn' 'varformat' 'vartype'          };
defaults = {100         sqrt(eps)   eps^(3/4)   []        []       ...
            []          []          'beta'      'hessian'          };
[eid emsg   MaxIter     TolParam    TolObj      Param0    Covar0   ...
            EstMethod   OutFun      VarFormat   VarType        ] = ...
                statgetargs(okargs,defaults,varargin{:});
if ~isempty(eid)
    error(sprintf('stats:mvregress:%s',eid),emsg);
end

if ~isscalar(MaxIter)
    error('stats:mvregress:BadMaxIter','MAXITER must be a scalar.');
elseif MaxIter < 1
    MaxIter = 1;
end

if ~isempty(EstMethod)
    if ~ischar(EstMethod) || size(EstMethod,1)~=1
        EstMethod = [];
    else
        EstMethod = strmatch(lower(EstMethod),{'cwls' 'ecm' 'mvn'});
    end
    if isempty(EstMethod)
        error('stats:mvregress:BadAlgorithm',...
              'ALGORITHM must be ''cwls'', ''ecm'', or ''mvn''.');
    end
end
if ~(isempty(EstMethod) || ismember(EstMethod,1:3))
    error('stats:mvregress:BadAlgorithm',...
          'ALGORITHM must be ''cwls'', ''ecm'', or ''mvn''.');
end

% Check variance format
okvals = {'beta' 'full'};
if ~ischar(VarFormat) || size(VarFormat,1)~=1
    VarFormat = [];
else
    VarFormat = strmatch(lower(VarFormat),okvals);
end
if isempty(VarFormat)
    error('stats:mvregress:BadVarFormat',...
          'VARFORMAT must be ''beta'' or ''full''.');
end
okvals{1} = 'paramonly';  % internal code for 'beta'
VarFormat = okvals{VarFormat};

% Check variance type
okvals = {'hessian' 'fisher'};
if ~ischar(VarType) || size(VarType,1)~=1
    VarType = [];
else
    VarType = strmatch(lower(VarType),okvals);
end
if isempty(VarType)
    error('stats:mvregress:BadVarType',...
          'VARTYPE must be ''hessian'' or ''fisher''.');
end
VarType = okvals{VarType};

% Check inputs, ignoring NaN rows for mvn method (3)
if isempty(EstMethod) && ~any(isnan(Data(:)))
    EstMethod = 3;   % use faster method for cases where others are equivalent
end
[NumSamples, NumSeries, NumParams, Data, Design, goodrows] = ...
          statcheckmvnr(Data, Design, Param0, Covar0, isequal(EstMethod,3));

celldesign = iscell(Design);
if celldesign && (numel(Design) == 1)
    SingleDesign = true;
else
    SingleDesign = false;
end
if ~celldesign && NumSeries>1
    error('stats:mvregress:BadDesign',...
          'X must be a cell array if Y has multiple columns.')
end

% Step 2 - observability and ignorability tests
Count = sum(all(isnan(Data),2));
if ((NumSamples - Count) * NumSeries) <= max(NumParams, (NumSeries * (NumSeries + 1))/2)
    if ((NumSamples - Count) * NumSeries) <= NumParams
        error('stats:mvregress:InsufficientData', ...
            'Insufficient data to estimate either full or least-squares models.');
    elseif isempty(EstMethod) || EstMethod ~= 1  % other than cwls method
        if isempty(EstMethod)
            EstMethod = 1; % select cwls
        else
            warning('stats:mvregress:TryingReducedModel', ...
                'Insufficient data for ''ecm'' algorithm.  Trying ''cwls''.');
            if EstMethod == 2  % ecm method
                EstMethod = 1; % switch to cwls
            else
                MaxIter = 1;   % or for mvn method, do just 1 iteration
            end
        end
    end
end
if isempty(EstMethod)
    EstMethod = 2;  % ecm method
end

% Step 3 - initialization
if isempty(Param0)
    Param = zeros(NumParams, 1);
else
    Param = Param0;
end

if isempty(Covar0)        % setup covariance-weighted least-squares
    CovarCWLS = eye(NumSeries, NumSeries);
    Covar = CovarCWLS;
else
    CovarCWLS = Covar0;
    Covar = Covar0;
end

VarParam = [];

[CholCovarCWLS, CholState] = chol(CovarCWLS);
if CholState > 0
    warning('stats:mvregress:NonPosDefCovar', ...
        'Input covariance matrix is not positive-definite. Will use an identity matrix.');
    CovarCWLS = eye(NumSeries, NumSeries);
    CholCovarCWLS = eye(NumSeries, NumSeries);
end

Resid = nan(NumSamples, NumSeries);

if ~SingleDesign
    if iscell(Design)
        Xdesign = reshape(permute(cat(3,Design{:}),[1 3 2]),[NumSeries*NumSamples,NumParams]);
    else
        Xdesign = Design;
    end
else
    Xdesign = repmat(Design{1},NumSamples,1);
end

if ~isempty(OutFun)
    str.Covar = CovarCWLS;
    str.iteration = 0;
    str.fval = [];
    if OutFun(Param,str,'init')
        disp('Iterations terminated prematurely by user.');
        Resid = [];
        return
    end
end

% Step 4 - main loop
nans = isnan(Data);
partialrows = find(any(nans,2))';  % these rows already removed for mvn method
Count = sum(~all(isnan(Data),2));
seps = sqrt(eps);
for Iter = 1:MaxIter
    Z = Data;
    if ~isempty(partialrows)  % always empty for mvn method
        if celldesign && SingleDesign
            Mean = Design{1} * Param;
        else
            Means = reshape(Xdesign*Param, [NumSeries, NumSamples]);
        end

        % Step 5 - parameter combined E and CM step
        WarnState = warning('off','MATLAB:nearlySingularMatrix');
        for i = partialrows
            if ~SingleDesign
                Mean = Means(:,i);
            end

            [mX, mY, CXX, CXY, CYY] = ecmpart(Data(i,:), Mean, Covar);

            if ~isempty(mY)
                P = isnan(Data(i,:));
                Q = ~P;
                Y = Data(i,Q)';

                Z(i,P) = mX + CXY * (CYY \ (Y - mY));
                Z(i,Q) = Y;
            end
        end
        warning(WarnState);
    end

    A = reshape(CholCovarCWLS' \ reshape(Xdesign,NumSeries,NumSamples*NumParams),[NumSeries*NumSamples,NumParams]);
    B = reshape(CholCovarCWLS' \ Z', NumSeries*NumSamples,1);
    Param = A \ B;

    % Step 6 - combined E and CM step to estimate covariance parameters
    if celldesign && SingleDesign
        Mean = Design{1} * Param;
        Means = [];
    else
        Means = reshape(Xdesign*Param, [NumSeries, NumSamples]);
    end

    Covar0 = Covar;
    Covar = zeros(NumSeries,NumSeries);

    if SingleDesign
        Resid = Data - repmat(Mean',size(Z,1),1);
    else
        Resid = Data - Means';
    end
    CovAdj = zeros(NumSeries, NumSeries);

    if ~isempty(partialrows)  % always empty for mvn method
        WarnState = warning('off','MATLAB:nearlySingularMatrix');
        Z = zeros(1,NumSeries);
        for i = partialrows
            if ~SingleDesign
                Mean = Means(:,i);
            end

            [mX, mY, CXX, CXY, CYY] = ecmpart(Data(i,:), Mean, Covar0);

            if ~isempty(mY)
                CovAdj(:) = 0;
                if isempty(mX)
                    Z = Data(i,:);
                else
                    P = isnan(Data(i,:));
                    Q = ~P;
                    Y = Data(i,Q)';

                    Z(P) = mX + CXY * (CYY \ (Y - mY));
                    Z(Q) = Y;
                    CovAdj(P,P) = CXX - CXY * inv(CYY) * CXY';
                end

                Resid(i,:) = Z - Mean';

                Covar = Covar + CovAdj;
            end
        end
        warning(WarnState);
    end

    Covar = (Covar + Resid'*Resid) / Count;

    % Step 7 - evaluate objective and test for convergence
    %          update chol(covar) for non-lsq method only
    if EstMethod == 1  % cwls method
        Objective = statecmobj(Design,Data,Param,CovarCWLS,Resid,CholCovarCWLS);
    else
        [CholCovarCWLS, CholState] = chol(Covar);
        if CholState>0
            if EstMethod == 3
                eid = 'stats:statmvnrobj:NonPosDefCov';
            else
                eid = 'stats:statecmobj:NonPosDefCov';
            end
            error(eid, 'Covariance is not positive-definite.');
        end
        if EstMethod == 3     % mvn method
            Objective = statmvnrobj(Data,Design,Param,Covar,Resid,CholCovarCWLS);
        elseif EstMethod == 1 % cwls method
            Objective = statecmobj(Design,Data,Param,CovarCWLS,Resid,CholCovarCWLS);
        else
            Objective = statecmobj(Design,Data,Param,Covar,Resid,CholCovarCWLS);
        end
    end

    if Iter > 1
        TestObj = Objective - Objective0;
        TestParam = norm(Param - Param0)/sqrt(NumParams);

        EpsObj = TolObj * (1 + abs(Objective));
        EpsParam = TolParam * (seps + norm(Param));

        if ((TestObj >= 0.0) && (TestObj < EpsObj)) && (TestParam < EpsParam)
            break
        end
    end

    Objective0 = Objective;
    Param0 = Param;
    
    if ~isempty(OutFun)
        str.Covar = Covar;
        str.iteration = Iter;
        str.fval = Objective;
        if OutFun(Param,str,'iter')
            disp('Iterations terminated prematurely by user.');
            return
        end
    end

    if (Iter == MaxIter) && (MaxIter > 1) && (TolObj > 0.0) && (TolParam > 0.0)
        warning('stats:mvregress:EarlyTermination', ...
            'Maximum iterations completed. Convergence criterion not satisfied.');
        break
    end
    
    if EstMethod == 3 && NumSeries == 1
        break
    end
end

if ~isempty(OutFun)
    str.Covar = Covar;
    str.iteration = Iter;
    str.fval = Objective;
    OutFun(Param,str,'done'); 
end

% Restore resids to proper size
n = size(Resid,1);   % sample size used in fit
if nargout>=3 && ~all(goodrows);
    Resid = resizeresids(Resid,goodrows);
end

% Compute parameter var/covar if requested
if nargout>=4
    Info = statecmmvnrfish(Data,Design,Covar,VarType,VarFormat);
    VarParam = inv(Info);
end

% -------------------------------------------------
function [mX, mY, CXX, CXY, CYY] = ecmpart(z, m, C)
%ECMPART Partitioning function for missing data algorithms
% Private routine to partition a mean vector m and covariance matrix C
% according to the pattern of NaNs (missing values) in an input vector z.

P = isnan(z);
Q = ~P;

mX = m(P);
mY = m(Q);

CXX = C(P,P);
CXY = C(P,Q);
CYY = C(Q,Q);

function outr = resizeresids(inr,goodrows)
outr = nan(length(goodrows),size(inr,2));
outr(goodrows,:) = inr;
