function redimscreen_widehor
% Redimscreen orients the figure according to screen size while leaving the menu bar 
% of MATLAB visible.

%[left bottom width height]

%a= get (0, 'Screensize'); set(gcf, 'Position', [a(3)*0.3 a(4)*(0.51) a(3)*0.68 a(4)*0.42])
a= get (0, 'Screensize'); set(gcf, 'Position', [a(3)*0.02 a(4)*(0.51) a(3)*0.97 a(4)*0.42])

%axis fill
