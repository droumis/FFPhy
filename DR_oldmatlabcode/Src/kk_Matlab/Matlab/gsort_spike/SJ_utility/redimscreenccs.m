function redimscreenccs
% Redimscreen orients the figure according to screen size while leaving the menu bar 
% of MATLAB visible.
a= get (0, 'Screensize'); set(gcf, 'Position', [a(3)*0.3 a(4)*0.035 a(3)*0.4 a(4)*0.4])
%axis fill