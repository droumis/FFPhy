function redimscreen_2x2subplots
% Redimscreen orients the figure according to screen size while leaving the menu bar 
% of MATLAB visible.

%% HORIZONTAL OBLONG FIGURE FOR SSS_COMPARE


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
set(gcf, 'Position', [170 100 1200 850]);
