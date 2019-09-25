classdef nominal < categorical
%NOMINAL Create a nominal array.
%   B = NOMINAL(A) creates a nominal array from A.  A is a numeric, logical,
%   character, or categorical array, or a cell array of strings. NOMINAL
%   creates levels of B from the sorted unique values in A, and creates
%   default labels for them.
%
%   B = NOMINAL(A,LABELS) creates a nominal array from A, labelling the levels
%   in B using LABELS.  LABELS is a character array or cell array of strings.
%   NOMINAL assigns the labels to levels in B in order according to the sorted
%   unique values in A.
%
%   B = NOMINAL(A,LABELS,LEVELS) creates a nominal array from A, with possible
%   levels defined by LEVELS.  LEVELS is a vector whose values can be compared
%   to those in A using the equality operator.  NOMINAL assigns labels to each
%   level from the corresponding elements of LABELS.  If A contains any values
%   not present in LEVELS, the levels of the corresponding elements of B are
%   undefined.  Pass in [] for LABELS to allow NOMINAL to create default labels.
%
%   B = NOMINAL(A,LABELS,[],EDGES) creates a nominal array by binning the
%   numeric array A, with bin edges given by the numeric vector EDGES.  The
%   uppermost bin includes values equal to the rightmost edge.  NOMINAL
%   assigns labels to each level in B from the corresponding elements of
%   LABELS.  EDGES must have one more element than LABELS.
%
%   By default, an element of B is undefined if the corresponding element of A
%   is NaN (when A is numeric), an empty string (when A is character), or
%   undefined (when A is categorical).  NOMINAL treats such elements as
%   "undefined" or "missing" and does not include entries for them among the
%   possible levels for B.  To create an explicit level for those elements
%   instead of treating them as undefined, you must use the LEVELS input, and
%   include NaN, the empty string, or an undefined element.
%
%   You may include duplicate labels in LABELS in order to merge multiple
%   values in A into a single level in B.
%
%   See also ORDINAL, HISTC.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:32:16 $

    methods
        function b = nominal(varargin) %nominal(a,labels,levels,edges)
        b = b@categorical(varargin{:});
        end % nominal constructor
    end
end
