function [cl, nc]= findDiffCells(sel1, sel2, show)
% find cells that occur in in sel1 but NOT in sel2

if nargin<3 | isempty(show); show= 0; end

%[c1,n1]= collectCellList(sel1);
%[c2,n2]= collectCellList(sel2);

c1= collectAnaList(sel1);
c2= collectAnaList(sel2);
n1= length(c1.rat);
n2= length(c2.rat);

pf= allPlaceFields(sel1);

cl.select1= sel1;
cl.select2= sel2;
nc= 0;
for i=1:n1

    j= find(strcmp(c1.rat{i}, c2.rat) &...
        c1.cellnum(i,1)==c2.cellnum(:,1)&...
        c1.cellnum(i,3)==c2.cellnum(:,3)&...
        c1.cellnum(i,4)==c2.cellnum(:,4) );
    if isempty(j)
        nc=nc+1;
        cl.rat{nc}= c1.rat{i};
        cl.cellnum(nc,:)= c1.cellnum(i,:);
        cl.day(nc)= c1.day(i);
        cl.traj(nc)= c1.traj(i);
        cl.pf(nc,:)= pf(i,:);
    end
end

if show
    for i=1:nc
        fprintf(1, 'i=%2.d: %s [%d %d %d %2.d], day= %2.d, traj %d, x=[%.1f %.1f]', ...
            i, cl.rat{i}, cl.cellnum(i,:), cl.day(i), cl.traj(i), cl.pf(i,:));
%        for j=1:length(cl.traj{nc})
%            fprintf(1, ' %d', cl.traj{i}(j));
%        end
        fprintf(1, '\n');
    end
end
