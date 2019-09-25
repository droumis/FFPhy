% KS1s2d: 2-dimensional Kolmogorov-Smirnov test of a single sample against a null distribution, 
%         optionally randomized.  The approximate significance level is satisfactory for N>20 or so
%         and pr<0.20.
%
%     Usage: [pr,Dmax,id] = ks1s2d(pts,'nullfn',{narg1},{iter},{'randfn'},{narg2],{doplot},...
%                                  {a1},{a2},{a3},{a4},{a5},{a6},{a7},{a8},{a9},{a10})
%
%         pts =       [n x 2] matrix of point coordinates for observed sample.
%         'nullfn' =  name (in single quotes) of function that returns the null probabilities
%                       of occurrence for the four quadrants, centered at p0=(x0,y0):
%                           quadprob = nullfn(pt0,...)
%                       where:
%                         pt0 =      [m x 2] matrix of center coordinates.
%                         ... =      optional additional arguments.
%                         quadprob = [m x 4] matrix of corresponding quadrant
%                                      probabilies in the standard sequence: 
%                                      upper right, upper left, lower left, lower right.
%         narg1 =     optional number of additional arguments to be passed to 'nullfn'
%                       [default = 0].
%         iter =      optional number of randomization iterations.  If iter==0 [the default],
%                       the approximate p-value is estimated as for Press et al. (1992).
%         'randfn' =  optional name (in single quotes, required if iter>0) of function 
%                       that returns a set of k random points from the null distribution:
%                           randpts = randfn(k,...)
%                       where:
%                         k =       scalar number of points to be returned.
%                         ... =     optional additional arguments.
%                         randpts = [k x 2] matrix of point coordinates.
%         narg2 =     optional number of additional arguments to be passed to 'randfn'
%                       [default = 0].
%         doplot =    optional boolean variable indicating, if true, that a plot of the 
%                       randomized null distribution of Dmax is to be produced [default = 0].
%         a1,... =    narg1+narg2 (<=10) optional additional arguments to be passed to the
%                       functions 'nullfn' and 'randfn'.
%         -----------------------------------------------------------------------------------
%         pr =        probability of observed test-statistic value under the null hypothesis.
%         Dmax =      observed test-statistic value.
%         id =        vector (length 2) of coordinates of the data point corresonding to the
%                       maximum value of the test statistic, Dmax.
%

% Press, WH, SA Teukolsky, WT Vetterling, BP Flannery. 1992. Numerical recipes in C: 
%   the art of scientific computing.  Cambridge University Press.

% RE Strauss, 5/28/03

function [pr,Dmax,id] = ks1s2d(pts,nullfn,narg1,iter,randfn,narg2,doplot,...
                               a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
  if (~nargin) help ks1s2d; return; end;
  
  if (nargin < 3)  narg1 = []; end;
  if (nargin < 4)  iter = []; end;
  if (nargin < 5)  randfn = []; end;
  if (nargin < 6)  narg2 = []; end;
  if (nargin < 7)  doplot = []; end;
  
  if (isempty(narg1))  narg1 = 0; end;
  if (isempty(narg2))  narg2 = 0; end;
  if (isempty(iter))   iter = 0; end;
  if (isempty(doplot)) doplot = 0; end;
  
  if (iter==0 & isempty(randfn))
    error('  KS1S2D: randfn must be supplied for iter > 0');
  end;
  
  [N,P] = size(pts);
  if (P~=2)
    error('  KS1S2D: points must be 2-dimensional.');
  end;
  
  nullfn_str = [nullfn,'(pts'];                     % Setup function call to nullfn
  if (narg1 > 0)
    for i = 1:narg1
      nullfn_str = [nullfn_str,sprintf(',a%d',a(i))];
    end;
  end;
  nullfn_str = [nullfn_str,');'];
  
  randfn_str = [randfn,'(N'];                       % Setup function call to randfn
  if (narg2 > 0)
    for i = 1:narg2
      randfn_str = [randfn_str,sprintf(',a%d',a(i+narg1))];
    end;
  end;
  randfn_str = [randfn_str,');'];

  
  quadprob_null = eval(nullfn_str);                 % Get null quadrant probabilities
  quadprob_obs = ks1s2df(pts);                      % Get observed quadrant proportions
  [Dmax,id] = max(abs(quadprob_obs-quadprob_null)); % Find Dmax and its position
  
  if (iter==0)                                      % Approximate probability
    r = corr(pts);
    num = Dmax*sqrt(N);
    den = 1 + sqrt(1-r*r)*(0.25 - 0.75/sqrt(N));
    pr = ksprob(num/den);
    
  else                                              % Randomized probability
    Dnull = zeros(iter,1);
    Dnull(1) = Dmax;
    for it = 2:iter
      pts = eval(randfn_str);
      quadprob_null = eval(nullfn_str);               % Get null quadrant probabilities
      quadprob_obs = ks1s2df(pts);                    % Get randomized quadrant proportions
      Dnull(it) = max(abs(quadprob_obs-quadprob_null)); % Find Dmax and stash
    end;
    
    pr = sum(Dnull>=Dmax);                            % Randomized p-value
    
    if (doplot)                                     % Plot of null distribution
      figure;
      histgram(Dnull);
      putxlab('Dmax');
    end;
  end;

  return;
  