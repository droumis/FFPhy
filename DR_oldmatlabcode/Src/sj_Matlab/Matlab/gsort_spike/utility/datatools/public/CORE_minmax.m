function [xmn,xmx,mni,mxi] = CORE_minmax_(x)
%CORE_MINMAX       Core computational routine for MINMAX.
%   [XMN,XMX,MNI,MXI] = CORE_MINMAX(X) returns scalars such that
%   XMN = X(MNI) = min(X) and XMX = X(MXI) = max(X).  Ties in indices are
%   broken in favor of the lowest magnitude index.
%
%   CONDITIONS
%   ----------
%   X must be a real vector of type DOUBLE.  An N-D array X is treated as X(:).
%   Infinite values are allowed.   NaN's are ignored.

[xmn,mni] = min(x(:));
[xmx,mxi] = max(x(:));


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array = randn(2e6,1);
% tic; [minval,maxval,minind,maxind] = CORE_minmax(array);  t(1) = toc;
% tic; [mnval,mnind] = min(array);   [mxval,mxind] = max(array); t(2) = toc;
% printf('CORE_minmax took %5.3f sec and equivalent native code took %5.3f sec.', t(1), t(2));
% if (~isequal([minval,maxval,minind,maxind], [mnval,maxval,mnind,mxind]))
%     printf('The two calls did not produce the same results.');
% end
