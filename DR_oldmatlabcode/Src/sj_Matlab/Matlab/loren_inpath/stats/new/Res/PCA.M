% PCA: Default principal component analysis of data matrix X, using:
%            covariance matrix
%          3 components
%       1000 bootstrap iterations
%        95% confidence intervals on loadings and percents variance
%

% RE Strauss, 6/16/93
%   9/21/99 - changed sequence call for PCACOV.

[loadings,percvar,scores,not_used,CI_loadings,CI_percvar] ...
   = pcacov(X,3,[],[],1000);
whos
