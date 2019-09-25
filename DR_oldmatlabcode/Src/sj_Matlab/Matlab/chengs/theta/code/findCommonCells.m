function [cl, nc]= findCommonCells(sel1, sel2, show)
% find cells that occur in both selectid's

if nargin<3 | isempty(show); show= 0; end

[c1,n1]= collectCellList(sel1);
[c2,n2]= collectCellList(sel2);

cl.select1= sel1;
cl.select2= sel2;
nc= 0;
for i=1:n1

    j= find(strcmp(c1.rat{i}, c2.rat)' &...
        c1.cellnum(i,1)==c2.cellnum(:,1)&...
        c1.cellnum(i,3)==c2.cellnum(:,3)&...
        c1.cellnum(i,4)==c2.cellnum(:,4) );
    if ~isempty(j)
        nc=nc+1;
        cl.rat{nc}= c1.rat{i};
        cl.cellnum(nc,:)= c1.cellnum(i,:);
        cl.day(nc)= c1.day(i);
    end
end

if show
    for i=1:nc
        fprintf(1, 'i=%2.d: %s [%d x %d %2.d], day= %2.d\n', ...
            i, cl.rat{i}, cl.cellnum(i,[1,3,4]), cl.day(i));
    end
end
