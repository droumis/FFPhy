% TODEGS: Converts radians to degrees.
%
%     Usage: degs = todegs(rads)
%
%         rads = matrix of angles in radians.
%         -------------------------------------------------
%         degs = corresponding matrix of angles in degrees.
%

% RE Strauss, 2/26/00

function degs = todegs(rads)
  degs = rads .* 180 ./ pi;
  return;
