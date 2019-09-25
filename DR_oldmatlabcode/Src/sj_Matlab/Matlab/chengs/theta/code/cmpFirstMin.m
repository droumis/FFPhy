%function cmpTheta

setRoot

id= 'pf9';

sets{1}.label= 'novel1';
sets{1}.armId= 'novelArm';
sets{1}.plotcol= [1 0 0];

sets{2}.label= 'fam1';
sets{2}.armId= 'famArm';
sets{2}.plotcol= [0 0 0];


nmin=5;

load([root '/work/mvel_60s_pf9']);
for is=1:length(sets)
    if isempty(sets{is}); continue; end
    armId= sets{is}.armId;
    selectid= [id '-' armId];

    load(sprintf('%s/data/tetrodes-%s', root, selectid));
    load(sprintf('%s/data/analist-%s', root, selectid));
    load(sprintf('%s/work/coefAbs_%s', root, sets{is}.label));
    load(sprintf('%s/work/psana-60s-%s-%s', root, selectid, armId));
%    load(sprintf('%s/work/ripple/data_%s_%s_day1', ...  root, armId, id));

    cl= collectCellList(selectid);
    day2cl= find(cl.day==1);
    cl2day= nan*ones(length(cl.day), 1);
    cl2day(day2cl)= [1:length(day2cl)]';
    pf= allPlaceFields(selectid);

%for ia=1:2
%    sp= [];
%    for j=2
%        N= sum(nRips{ia}{j}(iepoch{ia}{j},:), 2);
%        n= sum(nSpikes{ia}{j}, 2);
%        spikerate{ia}{j}= nan*ones(size(N));
%        ind= N>= nmin;
%        spikerate{ia}{j}(ind)= n(ind)./ N(ind);
%        sp= [sp; n(ind)./ N(ind)];
%    end
%    SPRATE{3-ia}= sp;
%end


    oldrat= ''; oldd= -1;
    nc= length(coefAbs.val);
    coef{is}= coefAbs.val;
    for ic=1:nc
        i= coefAbs.ind(ic);
        if analist.day(i)~= 1; error('should not get here'); end
        rat= analist.rat{i}; d= analist.cellnum(i,1); e= analist.cellnum(i,2); 
        t= analist.cellnum(i,3); c= analist.cellnum(i,4);
        itet= find(strcmp(tetrodes.rat, rat)' & d== tetrodes.num(:,1) & e== tetrodes.num(:,2) & t== tetrodes.num(:,3));

        psmax{is}(:,ic)= ana.max(itet);
        mv{is}(:,ic)= mvel.(armId).(rat){d}{e};
        if ~strcmp(rat, oldrat) | d~= oldd
            load(sprintf('%s/%s/data2/spikedata%.2d.mat', root, rat, d));
    %        load(sprintf('%s/%s/data2/behavdata%.2d.mat', root, rat, d));
            oldrat= rat; oldd= d;
        end

        sd= spikedata{d}{e}{t}{c};
        ind= find(pf(i,1)<sd.linpos & sd.linpos<pf(i,2) & sd.time<ana.tmax(itet, analist.traj(i)+1) & sd.traj== analist.traj(i));
    %    nspikes(:,ic)= sum(behavdata{d}{e}.traj(sd.index(ind))==analist.traj(i));
        nspikes{is}(:,ic)= length(ind);

        icell= find(strcmp(cl.rat, rat)' & d==cl.cellnum(:,1) &...
            e==cl.cellnum(:,2) & t==cl.cellnum(:,3) & c==cl.cellnum(:,4));
%        sp_rip{is}(:,ic)= SPRATE{is}(cl2day(icell));
    end
end

opt.fitline=0; 

figure
opt.xstr= 'theta power 0-60s'; opt.ystr= {'correlation coef.', 'mean over 0-60s'};
opt.outname= 'corr_V_peak_day1';
for is=1:length(sets)
    if isempty(sets{is}); continue; end
    opt.plotcol= sets{is}.plotcol;
    linrel(psmax{is}, coef{is}, opt);
    hold on
end
hold off

figure
opt.xstr= 'n spikes'; opt.ystr= 'correlation coef';
opt.outname= 'corr_V_nspikes_day1';
for is=1:length(sets)
    if isempty(sets{is}); continue; end
    opt.plotcol= sets{is}.plotcol;
    linrel(nspikes{is}, coef{is}, opt);
    hold on
end
hold off

figure
opt.xstr= 'mean vel (cm/s) over 0-60s'; opt.ystr= 'correlation coef';
opt.outname= 'corr_V_mvel_day1';
for is=1:length(sets)
    if isempty(sets{is}); continue; end
    opt.plotcol= sets{is}.plotcol;
    linrel(mv{is}, coef{is}, opt);
    hold on
end
hold off



%@@ tmp

figure
opt.xstr= 'mean vel (cm/s) over 0-60s'; opt.ystr= 'correlation coef';
opt.outname= 'corr_V_mvel_large_day1';
for is=1:length(sets)
    if isempty(sets{is}); continue; end
    opt.plotcol= sets{is}.plotcol;
    ind= find(mv{is}>9);
    lmv{is}= mv{is}(ind);
    lcoef{is}= coef{is}(ind);
    linrel(lmv{is}, lcoef{is}, opt);
    hold on
end
ranksum(lmv{1},lmv{2})
ranksum(lcoef{1},lcoef{2})
hold off

keyboard
