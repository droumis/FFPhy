function [bcodes,alabels] = matchlevels(a,b)
%MATCHLEVELS Utility for matching levels of nominal arrays.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:32:21 $

alabels = a.labels;
convert = zeros(1,length(b.labels)+1,class(a.codes));
[tf,convert(2:end)] = ismember(b.labels,alabels);

noMatch = find(~tf);
if ~isempty(noMatch)
    if length(alabels)+length(noMatch) > categorical.maxCode
        error('stats:nominal:matchlevels:MaxNumLevelsExceeded', ...
              'Too many categorical levels.');
    end
   convert(noMatch+1) = length(alabels) + (1:length(noMatch));
   alabels = [alabels b.labels(noMatch)];
end

bcodes = reshape(convert(b.codes+1), size(b.codes));
