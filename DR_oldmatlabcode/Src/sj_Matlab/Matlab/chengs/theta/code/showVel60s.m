function showVel60s
% Calculate average velocity on novel and familiar arm over first 60s
%
% Almost same as showVel.m but without plotting code. Difference are indicated
% by @@.

id= 'pf9';
armId= {'famArm', 'novelArm'};

timeid= 'minocc';
timelabel= 'occupancy (min)';
maxt= 1;

global mv behavdata spikedata lindistpos vel info task
setRoot

mv= cell(3,1); 
for d=1:3; 
    for ia=1:2
        vtmp{d}{ia}= cell(maxt,1);
    end
end
    selectid= [id '-novel'];
[tet, epochs]= collectTet(selectid);
ne= length(epochs.rat);

ix= zeros(3,1);
for ie=1:ne
    rat= epochs.rat{ie}; d= epochs.num(ie,1); e= epochs.num(ie,2);
    if e~=4 & e~= 6; continue; end

    data2dir= fullfile(root,rat,'data2');
    datadir= fullfile(root,rat,'data');
    loadVar(data2dir, 'info', 0);
    loadVar(datadir, 'task', 0, rat, 1);

    loadVar(data2dir, 'behavdata', d);
    bd= behavdata{d}{e};
%            loadVar(datadir, 'lindistpos', d, rat, 1);
%            bd.traj= lindistpos{d}{e}.estinfothetavel;
%            bd.traj= lindistpos{d}{e}.estinfo;
    bd.traj(bd.ripple==1)= -1;
    loadVar(datadir, 'vel', d, rat, 1);
    for ia=1:length(armId)
        vtmp= auxGetMeanVel(rat, d, e, bd, armId{ia}, timeid, maxt);
        mvel.(armId{ia}).(rat){d}{e}= vtmp(1); %@@ different from showVel.m!
    end
end

save(sprintf('%s/work/mvel_60s_%s.mat', root, id), 'mvel');


function vtmp= auxGetVel(rat, d, e, bd, selectid, timeid, maxt)
% Select times only in either the novel arm (novelArm) or familiar
% satellite arm (famArm) during the novel exposure epochs (4,6, ....),
% or all (novel).
% OR the satellite arms in the familiar configuration (fam).
%
%   selectid:    'fam', 'novel', 'novelArm' or 'famArm'

% define satellite arm
%pad= 5;
%armlength= 65;
%minspikes= 200;
minspikes= 1;


global task info vel

mv= nan*ones(maxt,1);
v= interp1(vel{d}{e}.data(:,1), vel{d}{e}.data(:,2), bd.time);
v= v*task{d}{e}.pixelsize;

switch selectid
case 'fam'
case 'novel'
case 'novelArm'
case 'famArm'
otherwise error('specify valid selectid as option');
end


% determine which arm or trajs are novel
% 1 is always home arm, 3 is familiar left satellite, 7 is familiar right
if task{d}{e}.task(2)== 3 & task{d}{e}.task(3)== 7
    error('should not get here: no novel arms in this epoch');
end

switch selectid
case 'novelArm'
    if task{d}{e}.task(2)~= 3; trajsel= [0 1]; else trajsel= [2 3]; end
case 'famArm'
    if task{d}{e}.task(3)~= 7; trajsel= [0 1]; else trajsel= [2 3]; end
case 'novel'
    trajsel= [0:3];
case 'fam'
    trajsel= [0:3];
end

%    minpos= info{d}{e}.centerlinpos+pad;
%    maxpos= minpos+armlength;
%    if(maxpos >= info{d}{e}.maxlinpos); error('inconsitency'); end

%@@
minpos= info{d}{e}.centerlinpos;
maxpos= info{d}{e}.maxlinpos;

for itraj=1:length(trajsel)
    traj= trajsel(itraj);
    inarm= find(minpos<=bd.linpos & bd.linpos<=maxpos & bd.traj==traj);
%    inarm= find(bd.linpos>=minpos & bd.traj==traj);
%    inarm= find(bd.traj==traj);
    t= getTimes(timeid, rat, d, e, traj, maxt);
    oldt= bd.time(1);
    for it=1:min(length(t),maxt)
        ind{itraj}{it}= inarm(oldt<=bd.time(inarm) & bd.time(inarm)<t(it));
%        disp(length(ind{itraj}{it}));
        oldt= t(it);
    end
end

for it=1:maxt
    vtmp{it}= v([ind{1}{it}; ind{2}{it}]);
end

function mv= auxGetMeanVel(rat, d, e, bd, selectid, timeid, maxt)

vtmp= auxGetVel(rat, d, e, bd, selectid, timeid, maxt);
for it=1:maxt
    % require at least 200ms of data to estimate mean velocity
    if length(vtmp{it})<100; continue; end
    mv(it)= nanmean(vtmp{it});
end

