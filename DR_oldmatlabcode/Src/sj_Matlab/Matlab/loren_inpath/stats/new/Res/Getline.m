% GETLINE: Given the coordinates of two points, finds the slope and intercept
%          of the line connecting the points.
%
%     Usage: [slope,intercept] = getline(P1,P2)
%
%         P1 = 2-element vector of [x,y] coordinates of 1st point.
%         P2 = 2-element vector of [x,y] coordinates of 2nd point.
%         ---------------------------------------------------------
%         slope =     slope of line connecting points.
%         intercept = intercept of line connecting point.
%

function [slope,intercept] = getline(P1,P2)
  X1 = P1(1);
  X2 = P2(1);
  Y1 = P1(2);
  Y2 = P2(2);

  slope = (Y2-Y1)./(X2-X1);
  intercept = Y1 - slope .* X1;

  return;
  