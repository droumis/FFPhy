% LINELABL: Plots a text string next to a pre-existing line, based on orientation of line.
%           Does not actually plot the line or endpoints of the line.  
%
%     Usage: linelabl(endpoints,label,position,max_x_kern,max_y_kern,{fontsize})
%
%         endpoints -   [2 x 2] matrix of X (col 1) and Y (col 2) coordinates of the two 
%                         endpoints of the line (rows).
%         label -       character string (row vector) to be printed.
%         position -    2-element vector indicating coordinates of lower-left corner of 
%                         label.
%         max_x_kern -  max horizontal displacement of label from position, in real axis 
%                         units.
%         max_y_kern -  max vertical displacement of label from position, in real axis 
%                         units.
%         fontsize -    optional font size for label [default = 10];
%

% RE Strauss

function linelabl(endpoints,label,position,max_x_kern,max_y_kern,fontsize)
  if (nargin < 6)
    fontsize = 10;
  end;

  p = endpoints(1,:);               % Find angle of line with respect to horizontal
  q = p + [1e6,0];
  r = endpoints(2,:);
  theta = angl(q,p,r);

  if (theta < 0)                    % Adjust angle to upper cartesian quadrants
    theta = theta + 2*pi;
  end;
  if (theta > pi)                   
    theta = theta - pi;
  end;

  xkern = max_x_kern * sin(theta);
  ykern = max_y_kern * cos(theta);

  hold on;
  h =text(position(1)+xkern,position(2)+ykern,label);
  set(h,'fontsize',fontsize);

  return;

