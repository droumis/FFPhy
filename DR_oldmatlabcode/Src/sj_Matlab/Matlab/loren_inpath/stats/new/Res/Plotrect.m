% PLOTRECT: Plots a single rectangle, given point coordinates of two opposing 
%           corners of the rectangle and its color.
%
%     Usage: plotrect(C1,C2,color)
%
%           C1 =    2-element vector of point crds of one corner.
%           C2 =    2-element vector of point crds of opposite corner.
%           color = optional color/linetype string [default = 'k'].
%

% RE Strauss, 6/12/99
%   9/3/99 - changed plot color for Matlab v5.

function plotrect(C1,C2,color)
  if (nargin<3) color = 'k'; end;

  rect = [C1; C2(1) C1(2); C2; C1(1) C2(2); C1];
  plot(rect(:,1),rect(:,2),color);

  return;

