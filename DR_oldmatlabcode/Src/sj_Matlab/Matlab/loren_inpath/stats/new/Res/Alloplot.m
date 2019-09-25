% ALLOPLOT: Produces a separate horizontal bar-chart depiction of allometric 
%           coefficients for each of g groups.
%
%     Usage: alloplot(allom,{labl})
%
%           allom = [p x g] matrix of allometric coefficients for p characters
%                     and g groups.
%           labl =  optional [p x c] matrix of character labels, for p characters 
%                     of maximum label-length c.
%

% RE Strauss, 12/16/97
%   9/20/99 -   update handling of null input arguments.
%   11/25/99 -  changed plot colors for Matlab v5.

function alloplot(allom,labl)
  if (nargin < 2) labl = []; end;

  [p,g] = size(allom);

  if (~isempty(labl))                     % Flip labels upside-down
    labl = flipud(labl);
  end;

  for gi = 1:g                            % Cycle thru groups
    ax = flipud(allom(:,gi));
    ay = [1:p]';
%    [y,x] = bar(ay,ax);
    h = barh(ay,ax);
    x = get(h,'XData');
    y = get(h,'YData');
xy = [x y]


    ix = (x==0);
    x(ix) = ones(sum(ix),1);
    xpos = max(x) + 0.03*range(x);
    leny = length(y);
    y(1) = y(2) - 0.3*(y(4)-y(3));
    y(leny) = y(leny-1) + 0.3*(y(leny-2)-y(leny-3));
    dy = 0.05*(y(leny) - y(1));

    figure;
    plot(x,y,'k');
    hold on;
    plot([1 1],[y(1) y(leny)],'k');
    putbnd(x,y);
    set(gca,'Ytick',[]);                    % Suppress y-axis labels and tick marks
    set(gca,'Ycolor','w');                  % Make y-axes invisible
    for pi = 1:p
      if (isempty(labl))
        text(xpos,pi,int2str(p-pi+1));
      else
        text(xpos,pi,labl(pi,:));
      end;
    end;
    putxlab('Allometric Coefficients');
    hold off;
  end;


  return;

