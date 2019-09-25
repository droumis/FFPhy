function redimscreen_2horsubplots
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
%set(gcf, 'Position', [a(3)*0.15 a(4)*0.15 a(3)*0.43 a(4)*0.5])
set(gcf, 'Position', [144 429 1622 588])