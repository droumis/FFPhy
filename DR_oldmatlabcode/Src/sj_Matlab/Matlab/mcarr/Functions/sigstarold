function h = sigstar(x1,x2,y,yh,starsym)

h{1} = line([x1; x1; x2; x2],[y; y+yh; y+yh; y],'Color','k','LineWidth',2);

midx = mean([x1,x2]);

if nargin < 5
  starsym = '*';
end

h{2} = text(midx,y+yh,starsym,'FontSize',18);
set(h{2},'horizontalalignment','center','verticalalignment','bottom');
