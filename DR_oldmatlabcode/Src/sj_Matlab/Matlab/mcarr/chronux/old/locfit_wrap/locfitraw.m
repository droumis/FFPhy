function [x,y,e]=locfitraw(varargin)
% locfitraw locfit helper function to call from matlab
%  
%  Usage: [x,y,e]=locfitraw( data )  {most basic usage, all defaults}
%
%    Additional arguments are attached as name-value pairs, ie:
%    [x,y,e]=locfitraw( data, 'alpha',[0.7,1.5] , 'family','rate' , 'ev','grid' , 'mg',100 ); 
%
%====================================================================
%
%  Argument types:
%
%      The first set of arguments ('x', 'y', 'weights', 'cens', and
%      'base') specify the regression variables and associated
%      quantities.
% 
%      Another set ('scale', 'alpha', 'deg', 'kern', 'kt', 'acri' and
%      'basis') control the amount of smoothing: bandwidth, smoothing
%      weights and the local model.
% 
%      'deriv' and 'dc' relate to derivative (or local slope) estimation.
% 
%      'family' and 'link' specify the likelihood family.
% 
%      'xlim' and 'renorm' may be used in density estimation.
% 
%      'ev', 'flim', 'mg' and 'cut' control the set of evaluation points.
% 
%      'maxk',  'itype', 'mint', 'maxit' and 'debug' control the Locfit
%      algorithms, and will be rarely used.
% 
%      'geth' and 'sty' are used by other functions calling 'locfit.raw',
%      and should not be used directly.
% 
%=========================================================================
%
%  Arguments in detail:
% 
%        x: Vector (or matrix) of the independent variable(s). 
%      ******************************
%  NOTE:       The first argument is placed in the first function slot without a name...
%              All other arguments require 'name',value notation
%      ******************************
%
%        y: Response variable for regression models. For density
%           families, 'y' can be omitted. 
% 
%  weights: Prior weights for observations (reciprocal of variance, or
%           sample size). 
% 
%     cens: Censoring indicators for hazard rate or censored regression.
%           The coding is '1' (or 'TRUE') for a censored observation, and
%           '0' (or 'FALSE') for uncensored observations. 
% 
%     base: Baseline parameter estimate. If provided, the local
%           regression model is fitted as Y_i = b_i + m(x_i) + epsilon_i,
%           with Locfit estimating the m(x) term. For regression models,
%           this effectively subtracts b_i from Y_i. The advantage of the
%           'base' formulation is that it extends to likelihood
%           regression models. 
% 
%    scale: A scale to apply to each variable. This is especially
%           important for multivariate fitting, where variables may be
%           measured in non-comparable units. It is also used to specify
%           the frequency for 'ang' terms. If 'scale=F' (the default) no
%           scaling is performed. If 'scale=T', marginal standard
%           deviations are used. Alternatively, a numeric vector can
%           provide scales for the individual variables. 
% 
%    alpha: Smoothing parameter. A single number (e.g. 'alpha=0.7') is
%           interpreted as a nearest neighbor fraction. With two
%           componentes (e.g. 'alpha=c(0.7,1.2)'), the first component is
%           a nearest neighbor fraction, and the second component is a
%           fixed component. A third component is the penalty term in
%           locally adaptive smoothing. 
% 
%      deg: Degree of local polynomial. Default: 2 (local quadratic).
%           Degrees 0 to 3 are supported by almost all parts of the
%           Locfit code. Higher degrees may work in some cases. 
% 
%     kern: Weight function, default = '"tcub"'. Other choices are
%           '"rect"', '"trwt"', '"tria"', '"epan"', '"bisq"' and
%           '"gauss"'. Choices may be restricted when derivatives are
%           required; e.g. for confidence bands and some bandwidth
%           selectors. 
% 
%       kt: Kernel type, '"sph"' (default); '"prod"'. In multivariate
%           problems, '"prod"' uses a simplified product model which
%           speeds up computations. 
% 
%       acri: Criterion for adaptive bandwidth selection.
% 
%    basis: User-specified basis functions. See 'lfbas' for more details
%           on this argument.
% 
%    deriv: Derivative estimation. If 'deriv=1', the returned fit will be
%           estimating the derivative (or more correctly, an estimate of
%           the local slope). If 'deriv=c(1,1)' the second order
%           derivative is estimated. 'deriv=2' is for the partial
%           derivative, with respect to the second variable, in
%           multivariate settings.  
% 
%       dc: Derivative adjustment.  
% 
%   family: Local likelihood family; '"gaussian"'; '"binomial"';
%           '"poisson"'; '"gamma"' and '"geom"'. Density and rate
%           estimation families are '"dens"', '"rate"' and '"hazard"'
%           (hazard rate). If the family is preceded by a ''q'' (for
%           example, 'family="qbinomial"'), quasi-likelihood variance
%           estimates are used. Otherwise, the residual variance ('rv')
%           is fixed at 1. The default family is '"qgauss"' if a response
%           'y' is provided; '"density"' if no response is provided. 
% 
%     link: Link function for local likelihood fitting. Depending on the
%           family, choices may be '"ident"', '"log"', '"logit"',
%           '"inverse"', '"sqrt"' and '"arcsin"'. 
% 
%     xlim: For density estimation, Locfit allows the density to be
%           supported on a bounded interval (or rectangle, in more than
%           one dimension). The format should be 'c(ll,ul)' where 'll' is
%           a vector of the lower bounds and 'ur' the upper bounds.
%           Bounds such as [0,infty) are not supported, but can be
%           effectively implemented by specifying a very large upper
%           bound. 
% 
%   renorm: Local likelihood density estimates may not integrate exactly
%           to 1. If 'renorm=T', the integral will be estimated
%           numerically and the estimate rescaled. Presently this is
%           implemented only in one dimension. 
% 
%       ev: Evaluation Structure, default = '"tree"'. Also available are
%           '"phull"', '"data"', '"grid"', '"kdtree"', '"kdcenter"' and
%           '"crossval"'. 'ev="none"' gives no evaluation points,
%           effectively producing the global parametric fit. A vector or
%           matrix of evaluation points can also be provided. 
% 
%     flim: A vector of lower and upper bounds for the evaluation
%           structure, specified as 'c(ll,ur)'. This should not be
%           confused with 'xlim'. It defaults to the data range. 
%         
%       mg: For the '"grid"' evaluation structure, 'mg' specifies the
%           number of points on each margin. Default 10. Can be either a
%           single number or vector. 
% 
%      cut: Refinement parameter for adaptive partitions. Default 0.8;
%           smaller values result in more refined partitions. 
% 
%     maxk: Controls space assignment for evaluation structures. For the
%           adaptive evaluation structures, it is impossible to be sure
%           in advance how many vertices will be generated. If you get
%           warnings about `Insufficient vertex space', Locfit's default
%           assigment can be increased by increasing 'maxk'. The default
%           is 'maxk=100'. 
% 
%    itype: Integration type for density estimation. Available methods
%           include '"prod"', '"mult"' and '"mlin"'; and '"haz"' for
%           hazard rate estimation problems. The available integration
%           methods depend on model specification (e.g. dimension, degree
%           of fit). By default, the best available method is used. 
% 
%     mint: Points for numerical integration rules. Default 20. 
% 
%    maxit: Maximum iterations for local likelihood estimation. Default
%           20. 
% 
%    debug: If > 0; prints out some debugging information.
% 
%     geth: Don't use!  
% 
%      sty: Style for special terms ('left', 'ang' e.t.c.). Do not try to
%           set this directly; call 'locfit' instead. 
%    
%==========================================================================
%
% Requires windows since R-(D)COM is windows-specific
%  I am working on a platform-independent replacement
%
% Requires that Matlab-R link Matlab package be installed from
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=5051&objectType=file   
% file MATLAB_RLINK.zip
%
% Requires that R be installed see http://r-project.org first
% file rw1091.exe
%
% Requires that R locfit package be installed first
% From within R in menu do "Packages" then "Install from CRAN"
%
% Requires that R-(D)COM be installed first from  
% http://lib.stat.cmu.edu/R/CRAN/contrib/extra/dcom/
% (get latest EXE file approx 3 MB)
% file RSrv135.exe
%
% The above packages should come bundled with this software for convenience
% with the exception of locfit which is easiest to install from within R
%
% In values:
% 
%

% Check for toolboxes
if not(exist('putRdata','file'));
    fprintf('You need to install Matlab-R Link first (do: "help locfitraw" for info)\nThen Install R-(D)COM\nThen install R\nThen install locfit from within R\nOnly works on Windoze\n');
    return
end

%   
% Connect to R only if not done so already, never disconnect
global RCONNECTED;
if isempty( RCONNECTED )
  % Try the open command
  [status,msg] = openR;
  if status ~= 1
    disp(['Problem connecting to R: ' msg]);
    return
  end
  evalR('library("locfit")') % attach locfit library
  RCONNECTED = 1;
end


% Minimal input validation    
if nargin < 1
   error( 'At least one input argument required' );
end
if mod(nargin,2)==0
   error( 'Argument count must be odd' );
end

putRdata( 'xdata', varargin{1}(:) );
args = '';

n = 2;
while n < length(varargin)
    if isa(varargin{n+1},'char')
      args = sprintf( '%s,%s="%s"',args, varargin{n}, varargin{n+1} );
    else
      putRdata( sprintf('%sval',varargin{n}), varargin{n+1} );
      args = sprintf( '%s,%s=%sval',args, varargin{n}, varargin{n} );
    end
    n=n+2;
end

command=sprintf( 'fit<-locfit.raw( xdata %s )', args );
evalR( command );
evalR( 'out<-knots(fit,what=c("x","coef","nlx"))' );
%evalR( 'plot(fit)' );

out = getRdata( 'out' );
[x,ind]=sort(out(:,1),1);
y=out(ind,2);
e=out(ind,3);
%aic=getRdata('-2*$fit$dp$lk+2*$fit$dp$df1');


return;
