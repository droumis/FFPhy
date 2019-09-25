function redimscreen100(window)
% Redimscreen orients the figure according to screen size while leaving the menu bar 
% of MATLAB visible.

%[left bottom width height]

a= get (0, 'Screensize'); 

switch window
    
    case 100
        factor = 0.25;
    case 101
        factor = 0.5;
    otherwise
        factor = 0.8;
end

set(gcf, 'Position', [a(3)*factor a(4)*0.035 a(3)*0.5 a(4)*0.90])
