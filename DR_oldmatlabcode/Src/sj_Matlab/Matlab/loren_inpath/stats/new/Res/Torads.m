% TORADS: Converts degrees to radians.
%
%     Usage: rads = torads(degs)
%
%         degs = matrix of angles in degrees.
%         -------------------------------------------------
%         rads = corresponding matrix of angles in radians.
%

% RE Strauss, 2/26/00

function rads = torads(degs)
  rads = degs .* pi ./ 180;
  return;
