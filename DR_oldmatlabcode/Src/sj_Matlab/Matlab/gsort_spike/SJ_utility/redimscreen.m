function redimscreen
% Redimscreen orients the figure according to screen size while leaving the menu bar 
% of MATLAB visible.
a= get (0, 'Screensize'); 
set(gcf, 'Position', [a(1)*0.95 a(2) a(3) a(4)*0.90])

