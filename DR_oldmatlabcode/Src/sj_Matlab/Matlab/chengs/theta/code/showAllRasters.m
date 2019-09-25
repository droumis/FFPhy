function showAllRasters(selectid, nsel)
%function showAllRasters(selectid, nsel)
%  selectid: which selection of cells to plot
%  nsel:     which placefields within selection to plot

maxt= 15;

timeid= 'passes';

setRoot;
load([root '/data/analist-' selectid])
nC= length(analist.rat);

novelDay= 0;
if nargin<2; 
    nsel=[1:nC]; 
    novelDay= 1;
end
Nsel= length(nsel);
sp= allSpikeData(selectid, 0, nsel);
pf= allPlaceFields(selectid);
tm= allTimes(timeid, selectid);

oldrat= ''; oldd= -1;
nx= ceil(sqrt(Nsel)); ny= ceil(Nsel/nx);
npass= 0;
nspikes= 0;

figure
for iC= 1:nC
    isel= find(iC== nsel);
    if isempty(isel); continue; end
    if novelDay>0 & analist.day(iC)~= novelDay; continue; end
    t= [0,tm{iC}];
    tsp= sp(iC).time;
    xsp= sp(iC).linpos;

%    subplot(nx, ny, isel);
    figure
    set(gcf, 'Name', ['raster-' num2str(iC)]);
    hold on 
    nt= min(length(t), maxt);
%    padding= diff(pf(iC,:))/10;
    padding= 0;
    xmin= pf(iC,1)-padding;
    xmax= pf(iC,2)+padding;
    for it=2:nt
        ind= find(t(it-1) < tsp & tsp<t(it) & ...
            xmin<xsp & xsp<xmax );
        x= xsp(ind); 
        n= length(x);
        npass= npass+1;
        nspikes(npass)= n;

        lh= line([x'; x'], [(it-1.4)*ones(1,n); (it-0.6)*ones(1,n)]);
        set(lh, 'Color', 'k');
        set(lh, 'LineWidth', 1);
        
    end
%    lh= line([xmin;xmax]*ones(1,nt), [1;1]*([1:nt]-0.5));
%    set(lh, 'Color', 'k');
%    lh= line([1;1]*pf(iC,:), [0;nt]*[1,1]);
%    set(lh, 'Color', 'k');
%    set(lh, 'LineStyle', ':');
    axis([xmin, xmax, 0, nt]);
    xlabel('position (cm)');
    ylabel('pass');
%    axis off

%    break
end
if length(nsel==1)
    print('-depsc2', ['rasters-' num2str(nsel)]);
else
    orient landscape
    print -depsc2 rasters
end

fprintf(1, 'totals: %d spikes, %d passes, %.2f spikes/pass, median= %d\n', ...
    sum(nspikes), npass, mean(nspikes), median(nspikes));
