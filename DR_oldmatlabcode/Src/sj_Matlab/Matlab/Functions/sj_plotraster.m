function h = sj_plotraster(x,y,h, clr, linewid)
% Shantanu - Removing figure no option - will plot on current figure
% Plot using "plot" instead of "line"

% [h] = plotraster(x,y,h,n, varargin)
%
% plots a vertical line of height h at each x, y coordinate in figure n.
% X & Y are Nx1 (column) vectors, h & n are scalars, suggest 0.8 for h.
%

%options: any plot properties like: 'color', 'k', 'linewidth', 4, 

if nargin < 3
    h = 0.8;
end

if nargin < 4
    clr = 'k';
end

if nargin < 5
    linewid = 2;
end

% Shantanu - use plot instead of line. Inefficient, but try for Illustrator purposes
for i = 1:length(x)   
    ao = y(i):0.1:y(i)+0.9;
    plot(x(i)*ones(size(ao)), ao, [clr '-'],'Linewidth',linewid, 'Markersize',5,'MarkerFaceColor','w');
end



% From "plotraster"
% ----------------
% if isempty(varargin)
%     options = {'color', 'k'};
%     options = {'Linewidth', 1};
% else
%     options = varargin;
% end
% 
% % if ~isempty(n)
% %     figure(n)
% % end
% 
% plotx = [x'; x'];
% ploty = [y'; y'+ h];
% h = line(plotx, ploty, options{1:end});
% 
% end
