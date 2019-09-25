function showCellHist(num)


if nargin < 1
    [d,e,t,c,ncells]= startCellList;
    while ~isempty(d)
        if e== 2
%        if d>1 & e== 2
            auxshow(d,e,t,c);
            pause
        end
        [d,e,t,c]= getNextCell;
    end
else
    for k= 1:size(num,1)
        d=num(k,1); e= num(k,2); t= num(k,3); c= num(k,4);
        auxshow(d,e,t,c);
    end
end

function auxshow(d,e,t,c);

fig=figure(1);
set(fig,'Position', [10 10 1010 310]);
global fmaux task tetloc behav occ trajocc
loadVar(fmaux.datadir, 'task', 0, fmaux.prefix, 1);
loadVar(fmaux.datadir, 'tetloc', 0, fmaux.prefix, 1);
loadVar(fmaux.datadir, 'behav', d, fmaux.prefix, 1);
loadVar(fmaux.datadir, 'occ', d, fmaux.prefix, 1);
loadVar(fmaux.datadir, 'trajocc', d, fmaux.prefix, 1);

%figure(fig);
clf
showbehav(behav, occ, trajocc, task, [d e t c]);
title(sprintf('%s %d %d %d %d, %s %s', fmaux.prefix, d, e, t, c, ...
	      tetloc{d}{t}.region, tetloc{d}{t}.layer));
