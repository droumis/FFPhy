function threeholepunch()
%
% THREEHOLEPUNCH
%
% Adjusts the printing size of the current figure to 
% accommodate holes for a three-ring binder, moving
% everything away from the left margin by a little.
% Should be called after calls to ORIENT.
%
% Copyright John Pezaris, 4/97, pz@caltech.edu.

units = get(gcf, 'PaperUnits');
set(gcf, 'PaperUnits', 'inches');
pos = get(gcf, 'PaperPosition');
orientation = get(gcf, 'PaperOrientation');

% The meaning of POS changes according to orientation,
% because of the implicit rotation.
if (strcmp(orientation, 'portrait'))
  dpos = 0.25 - pos(1);
  if (dpos > 0)
     pos(1) = pos(1) + dpos;
     pos(3) = pos(3) - dpos;
  end
elseif (strcmp(orientation, 'landscape'))
  pos(4) = pos(4) - 0.5;
end

set(gcf, 'PaperPosition', pos);
set(gcf, 'PaperUnits', units);
