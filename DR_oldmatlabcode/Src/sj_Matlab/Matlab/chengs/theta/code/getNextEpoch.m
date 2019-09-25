function [d,e]= getNextEpoch

global fmaux select

n= fmaux.currentCell;
od= select.cellnum(n,1);
oe= select.cellnum(n,2);

while n < fmaux.nCells
    n= n+1;
    d= select.cellnum(n,1);
    e= select.cellnum(n,2);
    
    if d~= od | e~= oe
	fmaux.currentCell= n;
	return
    end
end

% only get here, if no next epoch was found
fmaux.currentCell= n;
d=[];
e=[];
t=[];
c=[];

