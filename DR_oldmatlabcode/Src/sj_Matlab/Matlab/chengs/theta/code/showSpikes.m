function showSpikes(selectid, nsel)

maxt= 16;
timeid= 'passes';

load(['/bach/theta/data/analist-' selectid])
nC= length(analist.rat);
if nargin<2; nsel=[1:nC]; end
tm= allTimes(timeid, selectid);
sp= allSpikeData(selectid, 1, nsel);
pf= allPlaceFields(selectid);

Nsel= length(nsel);

oldrat= ''; oldd= -1;
npass= 0;
nspikes= 0;

for iC= 1:nC
    isel= find(iC== nsel);
    if isempty(isel); continue; end
    t= [0,tm{iC}];
    tsp= sp(iC).time;
    xsp= sp(iC).linpos;
    psp= sp(iC).theta;

    rat= analist.rat{iC}; num= analist.cellnum(iC,:)
    d=num(1); e=num(2); tet=num(3); c=num(4);

    fh= figure;
%    set(fh, 'Position', [357 54 604 892]);
    set(fh, 'Name', sprintf('%s [%d %d %d %d], no. %d', rat, d, e, tet, c, iC));
    nt= min(length(t), maxt);
    padding= diff(pf(iC,:))/10;
    xmin= pf(iC,1)-padding;
    xmax= pf(iC,2)+padding;
    for it=2:nt
        subplot(ceil((nt-1)/3), 3, it-1);
        hold on
        ind= find(t(it-1) < tsp & tsp<t(it) & ...
            xmin<xsp & xsp<xmax );
        x= xsp(ind); 
        phase= psp(ind);
        n= length(x);
        npass= npass+1;
        nspikes(npass)= n;

%        plot([x;x], [phase;phase+2*pi], 'k.')
        plot([x], [phase], 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
%        title(sprintf('pass %d, t= %5.2f',it-1, t(it)));
        title(sprintf('pass %d',it-1));
%        lh= line([1;1]*pf(iC,:), [0;4*pi]*[1,1]);
%        axis([xmin, xmax, 0, 4*pi]);
        lh= line([1;1]*pf(iC,:), [0;2*pi]*[1,1]);
        axis([xmin, xmax, 0, 2*pi]);
        set(lh, 'Color', 'k');
        set(lh, 'LineStyle', ':');
        axis square
    end
    orient tall
    print('-depsc2', sprintf('spikeplot-%s-%.2d', selectid, iC));

%    break
end

fprintf(1, 'totals: %d spikes, %d passes, %.2f spikes/pass, median= %d\n', ...
    sum(nspikes), npass, mean(nspikes), median(nspikes));
