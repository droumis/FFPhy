function [cnts,labs] = summary(a)
%SUMMARY Summary of a categorical array.
%   SUMMARY(A) displays the number of elements in the categorical array A
%   equal to each of A's possible levels.  If A contains any undefined
%   elements, the output also includes the number of undefined elements.
%
%   C = SUMMARY(A) returns counts of the number of elements in the categorical
%   array A equal to each of A's possible levels.  If A is a matrix or N-D
%   array, C is a matrix or array with rows corresponding to the A's levels.
%   If A contains any undefined elements, C contains one more row than the
%   number of A's levels, with the number of undefined elements in C(END) (or
%   C(END,:)).
%
%   [C,L] = SUMMARY(A) also returns the list of categorical level labels
%   corresponding to the counts in C.
%
%   See also CATEGORICAL/ISLEVEL, CATEGORCAL/ISMEMBER, CATEGORICAL/LEVELCOUNTS.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:14 $

c = levelcounts(a,1);
if nargout ~= 1
    labs = getlabels(a); labs = labs(:);
end
nundefs = sum(isundefined(a),1);
if any(nundefs(:) > 0)
    c = [c; nundefs];
    if nargout ~= 1
        labs = [labs; categorical.undefLabel];
    end
end

if nargout < 1
    if ndims(c) > 2
        tile = size(c); tile(1:2) = 1;
        labs = repmat(labs,tile);
    end
    c = [labs num2cell(c)];
    if size(a,2) == 1, c = c'; end
    str = evalc('disp(c)');
    str = str(1:end-1); % remove trailing newline
    % First find brackets containing numbers in any format, and preceeded by
    % whitespace -- those are the counts.  Replace those enclosing brackets
    % with spaces.  Then replace all quotes with spaces.
    str = regexprep(str,'(\s)\[([^\]]+)\]','$1 $2 ');
    str = regexprep(str,'''',' ');
    disp(str);
else
    cnts = c;
end
