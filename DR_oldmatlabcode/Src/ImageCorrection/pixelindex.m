function [ posind ] = pixelindex( pos, xdim, ydim )
% assigns an index to position data for tiff importers
%   Assumes tiff reader reads from the bottom left corner and proceeds row
%   by row until the upper right corner

posind=[pos(:,1)+pos(:,2).*xdim-xdim;];

% 
% 
% for i=1:size(pos,1);
%     currind=pos(i,1)+pos(i,2)*xdim-xdim;
%     posind=[posind;currind];
%     i=i+1;
% end

end

