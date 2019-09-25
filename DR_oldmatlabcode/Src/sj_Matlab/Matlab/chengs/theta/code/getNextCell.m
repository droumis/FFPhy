function [d,e,t,c]= getNextCell

global fmaux select

n= fmaux.currentCell;
if n < fmaux.nCells
	n= n+1;
    fmaux.currentCell= n;
    d= select.cellnum(n,1);
    e= select.cellnum(n,2);
    t= select.cellnum(n,3);
    c= select.cellnum(n,4);
    fprintf(1,'cell [%d %d %d %d], number %d\n',d,e,t,c,fmaux.currentCell);
else
    d=[];
    e=[];
    t=[];
    c=[];
end


