function matchaxis(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9)
% matchaxis - match limits accross multiple plots
%
%	MATCHAXIS sets the limits on all plots in the current figure to
%	the same values.  The new bottom (top) of each axis is the
%	minimum (maximum) of the bottoms (tops) of the existing axes.
%
%	MATCHAXIS(AXS), where AXS is a vector of axis handles, adjusts
%	the limits of all of the specified axes so that they are the same. 
%
%	MATCHAXIS(A1, A2 ...) is the same as MATCHAXIS([A1, A2 ...]);
%
%	See also AXES, AXIS, XAXIS, YAXIS.

% Copyright 1997, Maneesh Sahani maneesh@caltech.edu

ax = [];

if (nargin > 0) ax = [ax a0]; end
if (nargin > 1) ax = [ax a1]; end
if (nargin > 2) ax = [ax a2]; end
if (nargin > 3) ax = [ax a3]; end
if (nargin > 4) ax = [ax a4]; end
if (nargin > 5) ax = [ax a5]; end
if (nargin > 6) ax = [ax a6]; end
if (nargin > 7) ax = [ax a7]; end
if (nargin > 8) ax = [ax a8]; end
if (nargin > 9) ax = [ax a9]; end

if (isempty(ax)) ax = findobj(gcf, 'type', 'axes'); end

xlims = zeros(length(ax), 2);
ylims = zeros(length(ax), 2);
zlims = zeros(length(ax), 2);

for (i = 1:length(ax))
  xlims(i, :) = get(ax(i), 'Xlim');
  ylims(i, :) = get(ax(i), 'Ylim');
  zlims(i, :) = get(ax(i), 'Zlim');
end

Xlim = [min(xlims(:, 1)), max(xlims(:, 2))];
Ylim = [min(ylims(:, 1)), max(ylims(:, 2))];
Zlim = [min(zlims(:, 1)), max(zlims(:, 2))];

for (i = 1:length(ax))
  set (ax(i), 'Xlim', Xlim);
  set (ax(i), 'Ylim', Ylim);
  set (ax(i), 'Zlim', Zlim);
end
