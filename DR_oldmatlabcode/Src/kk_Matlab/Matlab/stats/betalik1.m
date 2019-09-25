function logL = betalik1(params,ld,l1d,n)
%BETALIK1 Obsolete function
%
%   Use BETALIKE in place of BETALIK1.

%Old help text follows.
%
%%BETALIK1 is a helper function. It is the same as BETALIKE, except it
% directly computes the betapdf without calling BETAPDF, also, saved some
% error checking and size checking.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.5.2.2 $  $Date: 2006/12/15 19:29:54 $

wmsgid = sprintf('stats:%s:Obsolete',mfilename);
warning(wmsgid, ...
    'The BETALIK1 function is now obsolete, use BETALIKE instead.\n');


a = params(1);
b = params(2);

logL = n*betaln(a,b)+(1-a)*ld+(1-b)*l1d;


