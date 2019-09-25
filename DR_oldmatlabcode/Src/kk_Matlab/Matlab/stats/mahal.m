function d = mahal(Y,X);
%MAHAL Mahalanobis distance.
%   D2 = MAHAL(Y,X) gives the Mahalanobis distance (in squared
%   units) of each point in Y from the sample in X. 
%   
%   X and Y must have the same number of columns, but can have
%   different numbers of rows.  X must have more rows than columns.
%
%   Example:  Generate some highly correlated X data.  The Y values
%             with equal coordinate values are much closer to X as
%             defined by Mahalanobis distance compared to the ones
%             with opposite coordinate values, even though they are
%             all equidistant from the mean using Euclidean distance.
%
%      x = mvnrnd([0;0], [1 .9;.9 1], 100);
%      y = [1 1;1 -1;-1 1;-1 -1];
%      mahal(y,x)

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.11.4.2 $  $Date: 2004/07/28 04:38:36 $

[rx,cx] = size(X);
[ry,cy] = size(Y);

if cx ~= cy
   error('stats:mahal:InputSizeMismatch',...
         'Requires the inputs to have the same number of columns.');
end

if rx < cx
   error('stats:mahal:TooFewRows',...
         'The number of rows of X must exceed the number of columns.');
end
if any(imag(X(:))) | any(imag(Y(:)))
   error('stats:mahal:NoComplex','MAHAL does not accept complex inputs.');
end

m = mean(X,1);
M = m(ones(ry,1),:);
C = X - m(ones(rx,1),:);
[Q,R] = qr(C,0);

ri = R'\(Y-M)';
d = sum(ri.*ri,1)'*(rx-1);

