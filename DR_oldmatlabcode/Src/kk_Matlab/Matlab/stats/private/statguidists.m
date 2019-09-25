function dists = statguidists
%STATGUIDISTS Distribution data for RANDTOOL and DISTTOOL

% The name field is used both to display the distribution, and as an argument
% to the CDF, PDF, ICDF, and RANDOM functions.  The rvname field is used
% internally to identify a distribution given its position in this structure,
% and is also used in RANDTOOL as a workspace variable name.
%
% The pmin and pmax fields give the theoretical limits for the parameters,
% while the plo(>=pmin) and phi(<=pmax) fields give the currently enforced
% limits.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/11/18 14:28:51 $

i = 0;
i = i+1; dists(i) = struct('name','Beta', ...
                           'rvname','betarv', ...
                           'discrete', 0, ...
                           'paramnames',{{'A','B'}}, ...
                           'parameters', [2 2], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [0 0], ...
                           'phi', [4 4], ...
                           'plo', [0.5 0.5], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Binomial', ...
                           'rvname','binorv', ...
                           'discrete', 1, ...
                           'paramnames',{{'Trials','Prob'}}, ...
                           'parameters', [10 0.5],  ...
                           'pmax', [Inf 1], ...
                           'pmin', [0 0], ...
                           'phi', [10 1], ...
                           'plo', [1 0], ...
                           'intparam', [1 0]);
i = i+1; dists(i) = struct('name','Chisquare', ...
                           'rvname','chi2rv', ...
                           'discrete', 0, ...
                           'paramnames',{{'df'}}, ...
                           'parameters', 2, ...
                           'pmax', Inf, ...
                           'pmin', 0, ...
                           'phi', 10, ...
                           'plo', 1, ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Discrete Uniform', ...
                           'rvname','unidrv', ...
                           'discrete', 1, ...
                           'paramnames',{{'Number'}}, ...
                           'parameters', 20, ...
                           'pmax', Inf, ...
                           'pmin', 0, ...
                           'phi', 20, ...
                           'plo', 1, ...
                           'intparam', [1 0]);
i = i+1; dists(i) = struct('name','Exponential', ...
                           'rvname','exprv', ...
                           'discrete', 0, ...
                           'paramnames',{{'Mu'}}, ...
                           'parameters', 1,  ...
                           'pmax', Inf, ...
                           'pmin', 0, ...
                           'phi', 2, ...
                           'plo', 0.5, ...
                           'intparam', [0]);
i = i+1; dists(i) = struct('name','Extreme Value', ...
                           'rvname','evrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'Mu','Sigma'}}, ...
                           'parameters', [0 1], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [-Inf 0], ...
                           'phi', [5 2], ...
                           'plo', [-5 .5], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','F', ...
                           'rvname','frv', ...
                           'discrete', 0, ...
                           'paramnames',{{'df1','df2'}}, ...
                           'parameters', [5 5], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [0 0], ...
                           'phi', [10 10], ...
                           'plo', [5 5], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Gamma', ...
                           'rvname','gamrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'A','B'}}, ...
                           'parameters', [2 2], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [0 0], ...
                           'phi', [5 5], ...
                           'plo', [2 2], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Generalized Extreme Value', ...
                           'rvname','gevrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'K','Sigma','Mu'}}, ...
                           'parameters', [0 1 0], ...
                           'pmax', [Inf Inf Inf], ...
                           'pmin', [-Inf 0 -Inf], ...
                           'phi', [.25 2 2], ...
                           'plo', [-0.5 0.5 -2], ...
                           'intparam', [0 0 0]);
i = i+1; dists(i) = struct('name','Generalized Pareto', ...
                           'rvname','gprv', ...
                           'discrete', 0, ...
                           'paramnames',{{'K','Sigma','Theta'}}, ...
                           'parameters', [0 1 0], ...
                           'pmax', [Inf Inf Inf], ...
                           'pmin', [-Inf 0 -Inf], ...
                           'phi', [.25 2 2], ...
                           'plo', [-0.25 0.5 -2], ...
                           'intparam', [0 0 0]);
i = i+1; dists(i) = struct('name','Geometric', ...
                           'rvname','georv', ...
                           'discrete', 1, ...
                           'paramnames',{{'Prob'}}, ...
                           'parameters', 0.5, ...
                           'pmax', 1, ...
                           'pmin', 0, ...
                           'phi', 0.99, ...
                           'plo', 0.25, ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Lognormal', ...
                           'rvname','lognrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'Mu','Sigma'}}, ...
                           'parameters', [0.5 0.25], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [-Inf 0], ...
                           'phi', [1 0.3], ...
                           'plo', [0 0.1], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Negative Binomial', ...
                           'rvname','nbinrv', ...
                           'discrete', 1, ...
                           'paramnames',{{'R','Prob'}}, ...
                           'parameters', [2 0.6], ...
                           'pmax', [Inf 1], ...
                           'pmin', [0 0], ...
                           'phi', [3 0.99], ...
                           'plo', [1 0.5], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Noncentral F', ...
                           'rvname','ncfrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'df1','df2','Delta'}}, ...
                           'parameters', [5 5 1], ...
                           'pmax', [Inf Inf Inf], ...
                           'pmin', [0 0 0], ...
                           'phi', [10 10 5], ...
                           'plo', [5 5 0], ...
                           'intparam', [0 0 0]);
i = i+1; dists(i) = struct('name','Noncentral T', ...
                           'rvname', 'nctrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'df','delta'}}, ...
                           'parameters', [5 1], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [0 -Inf], ...
                           'phi', [10 5], ...
                           'plo', [2 0], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Noncentral Chi-square', ...
                           'rvname', 'ncx2rv', ...
                           'discrete', 0, ...
                           'paramnames',{{'df','delta'}}, ...
                           'parameters', [5 1], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [0 0], ...
                           'phi', [10 5], ...
                           'plo', [2 0], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Normal', ...
                           'rvname', 'normrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'Mu','Sigma'}}, ...
                           'parameters', [0 1], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [-Inf 0], ...
                           'phi', [2 2], ...
                           'plo', [-2 0.5], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Poisson', ...
                           'rvname', 'poissrv', ...
                           'discrete', 1, ...
                           'paramnames',{{'Lambda'}}, ...
                           'parameters', 5, ...
                           'pmax', Inf, ...
                           'pmin', 0, ...
                           'phi', 5, ...
                           'plo', 1, ...
                           'intparam', [0]);
i = i+1; dists(i) = struct('name','Rayleigh', ...
                           'rvname', 'raylrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'B'}}, ...
                           'parameters', 2, ...
                           'pmax', Inf, ...
                           'pmin', 0, ...
                           'phi', 5, ...
                           'plo', 1, ...
                           'intparam', [0]);
i = i+1; dists(i) = struct('name','T', ...
                           'rvname', 'trv', ...
                           'discrete', 0, ...
                           'paramnames',{{'df'}}, ...
                           'parameters', 5, ...
                           'pmax', Inf, ...
                           'pmin', 0, ...
                           'phi', 10, ...
                           'plo', 2, ...
                           'intparam', [0]);
i = i+1; dists(i) = struct('name','Uniform', ...
                           'rvname', 'unifrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'Min','Max'}}, ...
                           'parameters', [0 1], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [-Inf -Inf], ...
                           'phi', [0.5 2], ...
                           'plo', [0 1], ...
                           'intparam', [0 0]);
i = i+1; dists(i) = struct('name','Weibull', ...
                           'rvname', 'weibrv', ...
                           'discrete', 0, ...
                           'paramnames',{{'A','B'}}, ...
                           'parameters', [1 2], ...
                           'pmax', [Inf Inf], ...
                           'pmin', [0 0], ...
                           'phi', [3 3], ...
                           'plo', [.5 .5], ...
                           'intparam', [0 0]);
