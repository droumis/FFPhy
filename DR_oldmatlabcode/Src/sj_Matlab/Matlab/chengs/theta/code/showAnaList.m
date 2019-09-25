function analist=showAnaList(arg1)

setRoot;
if ischar(arg1)
    selectid= arg1;
    load([root '/data/analist-' selectid]);
else
    analist= arg1;
end

pf= allPlaceFields(selectid);

na= length(analist.rat);
ja= 0;
%keyboard
for ia=1:na
    fprintf(1, '%3.d: %s [%2.d %d %d %2.d], day %d, traj= %d, x=[%.1f %.1f]\n', ...
        ia, analist.rat{ia}, analist.cellnum(ia,:), analist.day(ia),...
        analist.traj(ia), pf(ia,:));
end
