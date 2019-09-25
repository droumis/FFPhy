% nbars(x,y,iknots,prior,priorparam,burnin,sims,tau,c,fits,conf,bins,
%       verbose)
% 
% nbars.dll must be on the Matlab search path
% INPUT:
%  x: input data independent variable
%  y: input data dependent variable
%  iknots: initial number of knots e.g. 3
%  prior: prior form on number of knots must be 'Poisson', 'Uniform' or
%        'User'
%  priorparam: 'Uinform': 2D vector of range (min & max), e.g. [1 60].
%              'Poisson': scalar value of lambda, e.g. 6.0
%              'User': 2 column matrix of N rows (1st col = number of knots,
%                      2nd col = probability)
%  burnin: burnin iterations, e.g. 500
%  sims: sample iterations, e.g. 2000
%  tau: controls the spread for the knot proposal dist, e.g. 50.
%  c: reversible jump constant, e.g. 0.40
%  fits: number of grid points, e.g. 500
%  conf: confidence interval, e.g. 95
%  verbose: true or false indicate verbose output
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