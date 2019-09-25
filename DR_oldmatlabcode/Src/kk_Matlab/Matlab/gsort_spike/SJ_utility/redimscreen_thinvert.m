function redimscreen_thinvert(side)
% Redimscreen orients the figure according to screen size while leaving the menu bar 
% of MATLAB visible.

%% side: 0=left half of screen; 1=right half of screen

if nargin==0
    side=0;
end

a= get (0, 'Screensize'); 

switch side
    
    case 0
        factor = 0.02;
    case 1
        factor = 0.52;
    otherwise
        factor = 0.02;
end

%left bottom width height]
set(gcf, 'Position', [a(3)*factor a(4)*0.035 a(3)*0.25 a(4)*0.90])
