function h = weibplot(x)
%WEIBPLOT Obsolete function
%
%   Use WBLPLOT in place of WEIBPLOT.

%Old help text follows.
%
%WEIBPLOT Weibull probability plot.
%   H = WEIBPLOT(X) displays a Weibull probability plot of the  
%   data in X. For matrix, X, WEIBPLOT displays a plot for each column.
%   H is a handle to the plotted lines.
%   
%   The purpose of a Weibull probability plot is to graphically assess
%   whether the data in X could come from a Weibull distribution. If the
%   data are Weibull the plot will be linear. Other distribution types 
%   will introduce curvature in the plot.  

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.14.4.2 $  $Date: 2004/12/24 20:47:15 $

hndl = wblplot(x);
if nargout>0
   h = hndl;
end
