function showVel
% calculated average velocity on novel and familiar arm

rerun= 0;
onefig= 0;
calcid= 'moving2';
savefirstmin= 1;

id= 'pf9';
%id= 'all';
armId= {'famArm', 'novelArm'};

timeid= 'minocc';
timelabel= 'occupancy (min)';
maxt= 3;
%timeid= 'occ';
%timelabel= 'occupancy (s)';
%maxt= 300;

if savefirstmin
    timeid= 'minocc';
    timelabel= 'occupancy (min)';
    maxt= 1;
    rerun= 1;
end


global mv behavdata spikedata lindistpos vel info task
setRoot

fname= sprintf('%s/work/vel_%s_%s_%s.mat', root,  calcid,timeid,id);

if rerun
%    mv= {};
    mv= cell(3,1); 
    for d=1:3; 
        for ia=1:2
            vtmp{d}{ia}= cell(maxt,1);
        end
    end
%    selectid= [id '-novel'];
    selectid= id;
    [tet, epochs]= collectTet(selectid);
    ne= length(epochs.rat);

    ix= zeros(3,1);
    for ie=1:ne
        rat= epochs.rat{ie}; d= epochs.num(ie,1); e= epochs.num(ie,2);
        if e~=4 & e~= 6; continue; end

        novelday= epochs.day(ie);
        ix(novelday)= ix(novelday)+1;

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
            v= auxGetMeanVel(rat, d, e, bd, armId{ia}, timeid, maxt);
            if savefirstmin
                mvel.(armId{ia}).(rat){d}{e}= v(1);
            else
                for it=1:maxt
                    vtmp{novelday}{ia}{it}(ix(novelday),:)= v(it);
                end
            end

%            v= auxGetVel(rat, d, e, bd, armId{ia}, timeid, maxt);
%            for it=1:maxt
%                vtmp{novelday}{ia}{it}= [vtmp{novelday}{ia}{it}; v{it}];
%            end
        end
    end
    if savefirstmin
        save(sprintf('%s/work/mvel_60s_%s.mat', root, id), 'mvel');
        return;
    end
    save(fname, 'vtmp');
else
    load(fname);
end


    %auxPlotBars(mv, opt);
    plotcol= {[zeros(1,3); [1 0 0]]; [zeros(1,3); [0 1 0]]; [zeros(1,3); [0 0 1]]};
if onefig; 
    fh= figure; 
    set(fh,'position', [ 6   569   714   211]);
end
opt.xlabel= timelabel;
opt.ylabel='velocity (cm/s)';
opt.xticklabels= {'0-1', '1-2', '2-3'};
for d=1:3
%    if onefig; subplot(1,3,d); else figure; end
    if 1
        opt.title= sprintf('vel_%s_%s_%s_day%d', calcid, timeid,id,d);
        mv{d}=[] ;
        for ia=1:2
            for it=1:maxt
                vtmp{d}{ia}{it}= vtmp{d}{ia}{it}(vtmp{d}{ia}{it}>0); %@@
                mv{d}(it,ia)= nanmean(vtmp{d}{ia}{it});
%                dev{d}(it,ia)= nanstd(vtmp{d}{ia}{it});
                dev{d}(it,ia)= nanstd(vtmp{d}{ia}{it})/...
                    sqrt(sum(isfinite(vtmp{d}{ia}{it})));
            end
        end
%        hb= bar(mv{d});
        opt.plotcol= plotcol{d};
        auxBar(mv{d}, dev{d}, opt)
%        for ia=1:2; set(hb(ia), 'FaceColor', plotcol{d}(ia,:)); end
    else
        opt.plotcol= plotcol{d};
        auxPlotMeanDev([1:maxt], mv{d}, opt);
    end
    set(gca, 'ylim', [0 20]);
    if 0&~onefig
        xlabel(timelabel); ylabel('velocity (cm/s)');
        figname= sprintf('vel_%s_%s_%s_day%d', calcid, timeid,id,d);
        set(gcf, 'Name', figname);
        myprint([1.5 1], figname);
    end

    if 0
    sel{1}.label= 'novelArm';
    sel{2}.label= 'famArm';
    opt2.plot= 'cdf-p';
    opt2.outname= sprintf('%s_%s_%s_day%d', calcid,timeid,id, d);
    for it=1:maxt
        opt2.label= 'running speed (cm/s)';
        opt2.title= sprintf('vel_%dmin', it);
        auxCmpDist2({vtmp{d}{1}{it}, vtmp{d}{2}{it}}, sel, opt2);
    end
    end
end
if onefig
    subplot(1,3,1);
    ylabel('velocity (cm/s)');
    subplot(1,3,2);
    xlabel(timelabel);
    figname= sprintf('vel_%s_%s_%s', calcid,timeid,id);
    set(gcf, 'Name', figname);
    myprint([1.5 1], figname);
end


%v= {};
%for d=1:3
%    figure;
%    maxlen= 0;
%    for it=1:maxt
%        if length(vtmp{d}{ia}{it})>maxlen; maxlen= length(vtmp{d}{ia}{it}); end
%    end
%    for ia=1:2
%        v{ia}= nan*ones(500,maxlen);
%        for it=1:maxt
%            v{ia}(1:length(vtmp{d}{ia}{it}),it)= vtmp{d}{ia}{it};
%        end
%    end
%    auxPlotMeanDev([1:maxt], v, opt)
%    xlabel(timelabel); ylabel('velocity (cm/s)');
%    figname= sprintf('vel_%s_%s_%s_day%d', calcid, timeid,id,d);
%    set(gcf, 'Name', figname);
%    myprint([1.5 1], figname);
%end


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

if 1
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
else
    len= 60/0.002;
for itraj=1:length(trajsel)
    traj= trajsel(itraj);
    ontraj= find(minpos<=bd.linpos & bd.linpos<=maxpos & bd.traj==traj);
    ntraj= length(ontraj);
    lo= 1;
    for it=1:maxt; ind{itraj}{it}= []; end
    for it=1:maxt
        hi= lo+len-1;
        if hi>ntraj; 
            if ntraj>lo; hi= ntraj; else break; end
        end
        ind{itraj}{it}= ontraj(lo:hi);
        lo= hi+1;
    end
end
end

for it=1:maxt
    vtmp{it}= v([ind{1}{it}; ind{2}{it}]);
end

function mv= auxGetMeanVel(rat, d, e, bd, selectid, timeid, maxt)

vtmp= auxGetVel(rat, d, e, bd, selectid, timeid, maxt);
mv= nan*ones(maxt,1);
for it=1:maxt
    % require at least 200ms of data to estimate mean velocity
    if length(vtmp{it})<100; continue; end
    mv(it)= nanmean(vtmp{it});
end


function auxPlotBars(Zd, opt)
for d=1:3
    for j=1:2
        Zm(d,j)= nanmean(Zd{d,j});
        Zdev(d,j)= nanstd(Zd{d,j})/sqrt(sum(isfinite(Zd{d,j})));
%        if j>1
%             [p(d,j-1), h]= ranksum(Zd{d,j-1}, Zd{d,j});
%        end
    end
end
opt.xticklabels= {'day1'  'day2'  'day3'}; 
opt.plotcol={zeros(3,3), [[1 0 0]; [0 1 0]; [0 0 1]], [.7*ones(3,3)]};
auxBar(Zm, Zdev, opt, p);

