function [smooth,values,upper,lower]=rscb(in,w,nn,pts,lo,hi,pl)
%RSCB Runs locfit scb using R on a passed row array using the given parameters
% 
% RSCB( in, w, nn, pts, lo, hi )
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
%
%
% In values:
% 
% in: the input row array to smooth, contains a set of values 
% w: the smoothing width
% nn: nearest neighbor fraction for smoothing (0.1 is good)
% pts: number of evaluation points
% lo: lowest value to consider
% hi: highest value to consider
% pl: optional argument - make a plot if present and greater than zero
% 
% Out Values:
%
% smooth: a smoothed histogram of in, evaluated at the values points below
% values: the evaluation point values
% upper: the upper 95% confidence intervals for each point
% lower: the upper 95% confidence intervals for each point
%

% Minimal input validation    
if nargin < 6
   error( 'Not enough input arguments passed' );
end


%   
% Connect to R only if not done so already, never disconnect
global RCONNECTED;
if isempty( RCONNECTED )
  % Try the open command
  [status,msg] = openR;
  if status ~= 1
    disp(['Problem connecting to R: ' msg]);
  end
  evalR('library("locfit")') % attach locfit library
  RCONNECTED = 1;
end

% Put fitting parameters into R
putRdata( 'width', w' );
putRdata( 'nn', nn' );
putRdata( 'evalpts', pts );        
evalR( sprintf( 'flim=c(%f,%f)', lo, hi ) );
% bandwidth: variable and constant terms.  0.75*width/2 is how density.lf does it for a gaussian kernel
evalR( 'alpha=c(nn,0.75*width/2)' ); 

% send the data to R, transposed     
putRdata('data',in');


% set up for locfit.raw
evalR( 'data <- sort( data )' );
evalR( 'r <- range( data )' );        
%evalR( 'flim=c(r[1.]-width*0.75,r[2]+width*0.75) ); % alternative, data-derived limits - not based upon 'high' above

% run locfit.raw just like density.lf() and scb does
evalR( 'fit <- scb( data, ev="grid", mg=evalpts, flim=flim, alpha=alpha, kern="gauss", deg=0, link="ident", family="density", type=0 )' );
% pull evaluation points into matlab
evalR( 'values = fit$xev' );
values = getRdata( 'values' );

% pull smoothed data into matlab
evalR( 'smooth = fit$trans(fit$coef)' );
smooth = getRdata( 'smooth' );

% pull lower limits into matlab
evalR( 'lower = fit$trans(fit$lower)' );
lower = getRdata( 'lower' );

% pull upper limits into matlab
evalR( 'upper = fit$trans(fit$upper)' );
upper = getRdata( 'upper' );
    
if nargin > 6 && pl > 0
    % plot up results in matlab    
    clf % clear figure to prevent overplotting!
    plot( values, smooth, 'k-', values, [lower;upper], 'r:' );
    title( sprintf( 'Smoothed Results' ) );
    xlabel( 'time (ms)' );
    ylabel( 'Proportion' );
end
