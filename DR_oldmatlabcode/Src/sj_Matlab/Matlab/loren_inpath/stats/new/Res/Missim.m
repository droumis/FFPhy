% MISSIM: Simulate the effects of missing-data prediction methods (EM and PCA) 
%         for varying numbers of groups of observations and suites of variables.  
%         Creates a simulated data matrix and then calls the evalmiss() 
%         missing-data simulator.  Saves results to the file 'missim.txt'.
%
%   USAGE:  [goodresults,allbouts] = ...
%              missim(ngrp,nsuite,type,{pmiss},{iter},{ngrpmiss},{nsuitemiss},{wcorr},{bcorr})
%
%         ngrp =        vector of numbers of observations per group  
%                         (length = number of groups).
%         nsuite =      vector of numbers of characters per suite 
%                         (length = number of suites).
%         type =        estimation method: 0 for EM, 1 for PCA.
%         pmiss =       optional vector of proportions (0-1) or numbers (>1) of 
%                         missing data to be randomly introduced into data matrix
%                         [default = 0.01:0.01:0.50].
%         iter =        vector (length 3) containing optional numbers of 
%                         randomization iterations for each value of 'pmiss': 
%                         [total iterations, iterations per bout, maximum bouts]
%                         where number of bouts = (total iterations)/(iterations per 
%                         bout) but max number of bouts can be greater
%                         (default = [1000,200,10]).
%         ngrpmiss =    optional number of groups to contain missing values 
%                         [default = all].
%         nsuitemiss =  optional number of suites to contain missing values 
%                         [default = all].
%         wcorr =       optional within-suite correlation [default = 0.85].
%         bcorr =       optional between-suite correlation [default = 0.65].
%         -----------------------------------------------------------------------
%         goodresults = matrix containing results from bouts in which all 
%                         iterations were successful.  Also written to 
%                         'MissimGoodResults.txt' (see below).
%         allbouts =    matrix containing numbers of successful iterations
%                         per bout: [cur_pmiss, bout_number, percent_success]
%                         (see below).  Also written to 'MissimAllBouts.txt'.
%

% RE Strauss, 2/3/00
%   10/3/00 - add optional input values for within- and between-suite 
%               correlations.
%   10/5/00 - upgraded standardization and documentation.
%   3/6/01 -  move tofile to within 1:np loop.
%   11/5/02 - rewrite to accommodate changes in misseval().

function [goodresults,allbouts] = missim(ngrp,nsuite,type,pmiss,iter,ngrpmiss,nsuitemiss,wcorr,bcorr)

  if (nargin < 4) pmiss = []; end;
  if (nargin < 5) iter = []; end;
  if (nargin < 6) ngrpmiss = []; end;
  if (nargin < 7) nsuitemiss = []; end;
  if (nargin < 8) wcorr = []; end;
  if (nargin < 9) bcorr = []; end;

  if (isempty(pmiss))
    pmiss = [0.01:0.01:0.50];
  end;
  if (isempty(iter))
    iter = [1000,200,10];
  end;
  if (isempty(ngrpmiss))
    ngrpmiss = length(ngrp);
  end;
  if (isempty(nsuitemiss))
    nsuitemiss = length(nsuite);
  end;

  N = sum(ngrp);                        % Total number of observations
  P = sum(nsuite);                      % Total number of variables

  gr = [];
  ab = [];
  
  for it=1:length(pmiss)
    % Create simulated data matrix
    X = randsuit(ngrp,nsuite,[],[],wcorr,bcorr);   
    [h,mac,e] = homogen(X,1);

    % Evaluate performance of missing-data estimation method
    [goodresults,allbouts] = misseval(X,type,pmiss(it),iter,N,P);    
    
    o = ones(size(goodresults,1),1);
    goodresults = [it*o,pmiss(it)*o,h*o,mac*o,e*o,goodresults];
    gr = [gr; goodresults];
    tofile(goodresults,'MissimGoodresults.txt',-4,1);

    o = ones(size(allbouts,1),1);
    allbouts = [it*o,pmiss(it)*o,allbouts];
    ab = [ab; allbouts];
    tofile(allbouts,'MissimAllbouts.txt',-4,1);
  end;
  
  goodresults = gr;
  allbouts = ab;
    
  return;