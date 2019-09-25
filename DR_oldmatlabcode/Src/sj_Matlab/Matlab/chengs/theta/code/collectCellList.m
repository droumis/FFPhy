function [cl, nc]= collectCellList(selectid,show)
%function [cl, nc]= collectCellList(selectid)
%
%  Collect cell selections from all animals.
%

if nargin<2 | isempty(show); show= 0; end
rats= {'kyl', 'ter', 'fel', 'sta'};
nsel= length(rats);

cl= {}; nc=0;
setRoot;
for r=1:nsel
    load([root '/' rats{r} '/data2/select-' selectid])

    nr= size(select.cellnum,1);
    cl.cellnum(nc+1:nc+nr,:)=  select.cellnum;
    cl.day(nc+1:nc+nr,:)=  select.day;
    cl.newarm(nc+1:nc+nr,:)=  select.newarm;
    [cl.rat{nc+1:nc+nr}]= deal(rats{r});
    for i=1:nr
        for j=1:length(select.a{i})
            cl.traj{nc+i}(j)= select.a{i}{j}.traj;
            if isfield(select.a{i}{j}, 'linpos')
                cl.pf{nc+i}(:,j)= select.a{i}{j}.linpos;
            end
        end
    end
    nc= nc+nr;
end

if show
    for i=1:nc
        fprintf(1, 'i=%2.d: %s [%2.d %d %d %2.d], day= %2.d\n', ...
            i, cl.rat{i}, cl.cellnum(i,:), cl.day(i));
    end
end
