% pbars(x,y,logsp,iknots,prior,priorparam,burnin,sims,tau,c,fits,conf,
%       trials, bins, verbose)
%
% pbars.dll must be on the Matlab search path
% 
% INPUT:
%  X: input data independent variable (bin centers)
%  Y: input data dependent variable (spike counts, must be int32 type)
%  LOGSP: true or false use logspline to initialize knot placement
%  IKNOTS: initial number of knots e.g. 5
%  PRIOR: prior form on number of knots must be 'Poisson', 'Uniform' or
%        'User'
%  PRIORPARAM: 'Uinform': 2D vector of range (min & max), e.g. [1 15].
%              'Poisson': scalar value of lambda, e.g. 6.0
%              'User': 2 column matrix of N rows (1st col = number of knots,
%                      2nd col = probability)
%  BURNIN: burnin iterations, e.g. 100
%  SIMS: sample iterations, e.g. 1000
%  TAU: controls the spread for the knot proposal dist, e.g. 50.
%  C: reversible jump constant, e.g. 0.40
%  FITS: number of grid points, e.g. 150
%  CONF: confidence interval, e.g. 95
%  TRIALS: number of trials.
%  BINS: number of bins.
%  VERBOSE: true or false indicate verbose output
%
% OUTPUT:
% the output is written to a series of files
%  samp_mu: matrix length(x) by 'sims' fitted curves in columns
%  samp_mugrid: matrix 'fits' by 'sims' fitted curves in columns
%  samp_knots: under construction
%  samp_params: under construction
%  summ_mu: under construction
%  summ_mugrid: under construction
%  summ_params: under construction
%  summ_params: under construction