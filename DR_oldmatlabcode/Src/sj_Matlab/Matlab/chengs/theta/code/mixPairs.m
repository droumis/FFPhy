function [pl, np, nf]= mixPairs(sel1, sel2, show)
%function [pl, np, nf]= mixPairs(sel1, sel2)
% Find simultanously recorded pairs of cells (one cell from sel1, and one from
% sel2)

sel{1}= sel1;
sel{2}= sel2;
if nargin<3; show= 0; end

for i=1:2; [cl(i), nc(i)]= collectCellList(sel{i}); end

np= 0; nf= 0;
for ic=1:nc(1)
    rat= cl(1).rat{ic};
    num= cl(1).cellnum(ic,:); d=num(1); e=num(2); tet=num(3); c=num(4);
    traj= cl(1).traj{ic};

    % find simultaneous recorded cell
    ind= find(strcmp(rat, cl(2).rat') & ...
        (d== cl(2).cellnum(:,1)) & (e== cl(2).cellnum(:,2)) & ...
        (tet== cl(2).cellnum(:,3)) & (c== cl(2).cellnum(:,4)) );
    cl(2).cellnum(ind,:)= nan;
    ind= find(strcmp(rat, cl(2).rat') & ...
        (d== cl(2).cellnum(:,1)) & (e== cl(2).cellnum(:,2)) );
    if isempty(ind); continue; end

    n= length(ind);
    pl.cellnum(np+1:np+n, 1:4)= ones(n,1)*num;
    pl.cellnum(np+1:np+n, 5:8)= cl(2).cellnum(ind,:);
    [pl.rat{np+1:np+n}]= deal(rat);
    pl.day(np+1:np+n)= cl(2).day(ind);
    pl.newarm(np+1:np+n)= cl(2).newarm(ind);
    np= np+length(ind);
    nf= nan;
end

if show
    for i=1:np
        fprintf(1, '%s [%d %d %d %d], [%d %d %d %d]\n', ...
            pl.rat{i}, pl.cellnum(i,:));
    end
end
