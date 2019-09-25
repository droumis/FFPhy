function showISI(num, traj, nt, scale, vert)
%function showModelSnapshot(num, [traj, nt, scale, vert])
%
% scale=0:    use same colorscale for all plots

if nargin < 2; traj= [0 1 2 3]; end
if nargin < 3; nt= 4; end
if nargin < 4; scale= 0; end
if nargin < 5; vert= 'time'; end

if nargin < 1 | isempty(num)
    [d,e,t,c,ncells]= startCellList;
    while ~isempty(d)
        auxrun([d,e,t,c], traj, nt, scale, vert);
        [d,e,t,c]= getNextCell;
    end
else
    auxrun(num, traj, nt, scale, vert);
end

function auxrun(num, traj, nt, scale, vert)

global adaptest fmaux select
data=loadData(num);
loadVar('.','adaptest',num(1));
m= adaptest{num(1)}{num(2)}{num(3)}{num(4)}.model;

nt= 1;
[t, x,y]= auxEvalModel(m.isi, 0, nt, 500);

fh= figure;
for it=1:nt
    subplot(1,nt, it);
    semilogx(1000*x, y{1,it})

    % label plots
    ylabel('\lambda_T')
    xlabel('ISI (ms)'); 
end

myprint('large', 'isi');

function auxAppendTitle(str)
tstr= get(get(gca, 'title'), 'string');
if ~isempty(tstr); tstr= [tstr ', ']; end
title([tstr str]);
