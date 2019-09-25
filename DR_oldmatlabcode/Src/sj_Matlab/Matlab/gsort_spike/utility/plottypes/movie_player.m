function movie_player(mov, pauz)
%MOVIE_PLAYER      Quick and dirty way to watch an xyt matrix as a movie.
%   MOVIE_PLAYER(MOV) plays the pages of a 3-D matrix as a movie.  An
%   optional second argument is the interframe pause interval (sec).

if (nargin < 2)
    pauz = 0;
end

% use the 'pcolor' command if you want interpolated coloring, but if you
% don't 'image' is better because 'pcolor' won't show the last row/col.
% hand = pcolor(mov(:,:,1));
% shading interp
% %set(gcf,'Renderer','OpenGL');  % faster (but no colormap interpolation)
% set(hand,'CDataMapping','scaled');         % map to clim, . . . 

hand = imagesc(mov(:,:,1));
%set(gcf, 'Renderer', 'Painters', 'DoubleBuffer','on');

axis xy
set(gca,'CLim',[min(mov(:)) max(mov(:))]); % . . . set to the range of the data
%set(gca, 'CLim', [-3000, 3000])
h = title('');HGH
for k = 1:size(mov,3)
    set(hand,'CData',mov(:,:,k));
    set(h, 'String', ['Image #' num2str(k)]);
    pause(pauz);
    drawnow
end

