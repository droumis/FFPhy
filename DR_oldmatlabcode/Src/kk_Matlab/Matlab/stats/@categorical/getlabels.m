function labels = getlabels(a)
%GETLABELS Get level labels of a categorical array.
%   S = GETLABELS(A) returns the labels for the levels of the categorical
%   array A. S is a cell array of strings.
%
%   See also CATEGORICAL/SETLABELS.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:30:38 $

labels = a.labels;
