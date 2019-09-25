function h = plotraster_marker(x,y,h,n, varargin)
% [h] = plotraster(x,y,h,n, varargin)
%
% plots a vertical line of height h at each x, y coordinate in figure n.
% X & Y are Nx1 (column) vectors, h & n are scalars, suggest 0.8 for h.
%
%options: any plot properties like: 'color', 'k', 'linewidth', 4, 

if isempty(varargin)
    options = {'color', 'k'};
else
    options = varargin;
end
    
if ~isempty(n)
    figure(n)
end

plot(x,y,'bo','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','b');


% plotx = [x'; x'];
% ploty = [y'; y'+ h];
% h = line(plotx, ploty, options{1:end});

end
