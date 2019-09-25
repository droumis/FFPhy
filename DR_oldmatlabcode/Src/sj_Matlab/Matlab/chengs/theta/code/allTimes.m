function tm= allTimes3(timeid, selectid, maxt)
%function tm= allTimes3(timeid, selectid, maxt)
%
% Return times at which certain behavioral landmarks are reached for all cells.
% timeid: passes occ minocc acc minacc

load(['/home/chengs/theta/data/analist-' selectid])
tm= {};
nC= length(analist.rat);

for iC= 1:nC
    % set up data
    rat= analist.rat{iC};
    num= analist.cellnum(iC,:); d=num(1); e=num(2); tet=num(3); c=num(4);
    traj= analist.traj(iC);
    tm{iC}= getTimes(timeid, rat, d, e, traj, maxt);
end


