%function  showPerm(selectid, adaptid)
%function  showPerm(selectid, adaptid)

selectid= 'all-placefields2-HomeArm';
adaptid= 'adapt_xp015_t015';
permid= 'permstats-p05';

setRoot
load([root '/work/' permid '-' selectid '-' adaptid])
load([root '/data/analist-' selectid '-' adaptid])
%load(['collectstats-' selectid '-' adaptid])

vars= {}; opt= {};
rat= 'kyl';

nV=0;

nV=nV+1;
vars{nV}.name= 'MutualInfo';
vars{nV}.mean =1;
vars{nV}.pairdiff =1;

nV=nV+1;
vars{nV}.name= 'Integral2d';
vars{nV}.pairdiff=1;
vars{nV}.rel=1;
vars{nV}.mean =1;

nV=nV+1;
vars{nV}.name= 'MutualInfo';
vars{nV}.pair=1;
vars{nV}.mean =1;

nV=nV+1;
vars{nV}.name= 'Integral2d';
vars{nV}.pair=1;
vars{nV}.mean =1;

realvars= calcStats(realstats, vars, analist.iana);

for iV=1:nV
    realpro{iV}= realvars{iV};
    realpro{iV}.val= realvars{iV}.val(analist.proind,:);
    realpro{iV}.title=  'realpro';
    realretro{iV}= realvars{iV};
    realretro{iV}.val= realvars{iV}.val(analist.retroind,:);
    realretro{iV}.title=  'realretro';
end

cvopt=[];
%cvopt.title= 1;
cvopt.legend= 1;

%figure
%subplot(2,2,1); auxCmpVars(realpro{3}, realpro{4}, cvopt);
%cvopt.legend= 0;
%subplot(2,2,2); auxCmpVars(realretro{3}, realretro{4}, cvopt);

figure
subplot(2,2,1); auxCmpVars(realpro{1}, realpro{2}, cvopt);
cvopt.legend= 0;
subplot(2,2,2); auxCmpVars(realretro{1}, realretro{2}, cvopt);


nC= length(permstats);

for iV=1:nV
    ratpro{iV}= realpro{iV};
    ratpro{iV}.val= realvars{iV}.val(analist.([rat 'proind']),:);
    ratretro{iV}= realretro{iV};
    ratretro{iV}.val= realvars{iV}.val(analist.([rat 'retroind']),:);
end

permvars= {};
permpro= {};
permretro= {};
for iV=1:nV
    permvars{iV}= vars{iV};
    permvars{iV}.append =1;
    permvars{iV}.title = 'perm';

    permpro{iV}= vars{iV};
    permpro{iV}.append =1;
    permpro{iV}.title = 'permpro';
    ratpermpro{iV}= permpro{iV};

    permretro{iV}= vars{iV};
    permretro{iV}.append =1;
    permretro{iV}.title = 'permretro';
    ratpermretro{iV}= permretro{iV};
end

permind= ones(nC+1,1); % starting indices corresponds to permutations of iC
for iC= 1:nC
%    if iC==19; keyboard; end
    nP= length(permstats{iC});
    permind(iC+1)= permind(iC)+nP;
    iana= analist.iana(iC)*ones(nP,1);
    permvars= calcStats(permstats{iC}, permvars, iana);
    if find(analist.proind== iC)
        permpro= calcStats(permstats{iC}, permpro, iana);
        if find(analist.([rat 'ind' ])== iC) ratpermpro= calcStats(permstats{iC}, ratpermpro, iana); end
    elseif find(analist.retroind== iC)
        permretro= calcStats(permstats{iC}, permretro, iana);
        if find(analist.([rat 'ind' ])== iC) ratpermretro= calcStats(permstats{iC}, ratpermretro, iana); end
    else
        error('should not get here: traj has to be either pro or retro');
    end
end
%keyboard
%auxShowStats(permpro, opt, analist);
%auxShowStats(permretro, opt, analist);

subplot(2,2,3); auxCmpVars(permpro{1}, permpro{2}, cvopt);
subplot(2,2,4); auxCmpVars(permretro{1}, permretro{2}, cvopt);

figname= [root '/work/' permid '-' permpro{1}.name '-vs-' permpro{2}.name '-' selectid '-' adaptid];
print('-dpng', figname)
saveas(gcf, figname)

cdopt= [];
cdopt.nbins= 15;
%cdopt.title= 1;

% compare distributions
%permvars{1}.val vs realvars{1}.val
plotid= {'cdf', 'pdf', 'hist'};
for ip=1:length(plotid)
cdopt.plot= plotid{ip};
cdopt.legend= 1;
figure
subplot(2,2,1); auxCmpDist(realpro{2}, permpro{2}, cdopt);
cdopt.legend= 0;
subplot(2,2,2); auxCmpDist(realretro{2}, permretro{2}, cdopt);

subplot(2,2,3); auxCmpDist(realpro{1}, permpro{1}, cdopt);
subplot(2,2,4); auxCmpDist(realretro{1}, permretro{1}, cdopt);
figname= [root '/work/' permid '-' plotid{ip} '-' selectid '-' adaptid];
print('-dpng', figname)
saveas(gcf, figname)
end
%keyboard

%n= 1000;
%x= randn(n,1);
%realpro{1}.val= x;
%permpro{1}.val= x+.1;
%figure
%auxCmpDist(realpro{1}, permpro{1}, cdopt);


% plot for individual rat
%figure
%subplot(2,2,1); auxCmpDist(ratpro{2}, ratpermpro{2}, cdopt);
%cdopt.legend= 0;
%subplot(2,2,2); auxCmpDist(ratretro{2}, ratpermretro{2}, cdopt);

%subplot(2,2,3); auxCmpDist(ratpro{1}, ratpermpro{1}, cdopt);
%subplot(2,2,4); auxCmpDist(ratretro{1}, ratpermretro{1}, cdopt);

%figure
%iV= 2;
%plot([1,nC], [0 0], 'k-');  hold on
%for iC= 1:nC
%    plot(iC, realvars{iV}.val(iC), 'sk', 'markerfacecolor', [0 0 0]);
%    permval= permvars{iV}.val(permind(iC):permind(iC+1)-1);
%    plot(iC, permval, '.r');
%    line([iC, iC], [min(permval), max(permval)]);
%    if abs(min(permval)- max(permval)) > .3 ;        iC ;    end
%    if (realvars{iV}.val(iC) < min(permval) | ...
%    if abs(realvars{iV}.val(iC))>.3 & (realvars{iV}.val(iC) < min(permval) | ...
%        realvars{iV}.val(iC) > max(permval))
%        iC
%        plot(iC, realvars{iV}.val(iC), 'og');
%    end
%end

%keyboard
