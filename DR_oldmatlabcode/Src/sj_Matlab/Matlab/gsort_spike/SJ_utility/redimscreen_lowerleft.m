function redimscreen_lowerleft
% Redimscreen orients the figure according to screen size while leaving the menu bar 
% of MATLAB visible.

%[left bottom width height]

a= get (0, 'Screensize'); 
set(gcf, 'Position', [a(1)*0.95 a(2)*0.2 a(3)*0.5 a(4)*0.5])
