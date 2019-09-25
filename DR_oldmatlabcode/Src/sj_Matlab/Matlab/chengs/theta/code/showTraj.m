function showTraj(num)
%function showTraj(num)
%
% Show trajectory classifications.

global adaptest
if nargin < 1 | isempty(num)
    [d,e,t,c,ncells]= startCellList;
    while ~isempty(d)
        auxrun(d,e,t,c);
        [d,e,t,c]= getNextCell;
    end
else
    auxrun(num(1),num(2),num(3),num(4));
end

function auxrun(d,e,t,c)

global adaptest fmaux

data= loadData([d e t c]);
fh= figure;
figname= sprintf('%s [%d %d %d %d]', fmaux.prefix, d, e, t, c);
set(fh, 'Name', figname);

% mirror setup on x-axis
miny= min(data.ypos);
maxy= max(data.ypos);
data.ypos= (maxy+miny) - data.ypos;


col= [1 0 0; 0 1 0; 0 0 1; 1 1 0];
for itraj=1:4
    subplot(2,2,itraj);
    i= find(data.traj== itraj-1);
    plot(data.xpos, data.ypos, '.', 'MarkerEdgeColor', .7*[1 1 1])
    hold on
    plot(data.xpos(i), data.ypos(i), '.', 'MarkerEdgeColor', col(itraj,:));
    title(sprintf('traj: %d', itraj-1));
end

%keyboard

