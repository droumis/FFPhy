function redimscreen_portrait
% Redimscreen orients the figure according to screen size while leaving the menu bar 
% of MATLAB visible.

%% Normal Portrait Figure

if nargin==0
    side=0;
end

a= get (0, 'Screensize'); 

% switch side
%     
%     case 0
%         factor = 0.02;
%     case 1
%         factor = 0.52;
%     otherwise
%         factor = 0.02;
% end

%left bottom width height]
set(gcf, 'Position', [a(3)*0.05 a(4)*0.1 a(3)*0.35 a(4)*0.75])
