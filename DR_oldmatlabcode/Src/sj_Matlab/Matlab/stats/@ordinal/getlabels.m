function labels = getlabels(a)
%GETLABELS Get level labels of an ordinal array.
%   S = GETLABELS(A) returns the labels for the levels of the ordinal
%   array A. S is a cell array of strings.  S contains the labels ordered
%   according to the ordering of the levels of A.
%
%   See also ORDINAL/SETLABELS.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2006/12/15 19:32:26 $

labels = getlabels@categorical(a);
