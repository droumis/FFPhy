function showConcurrent(selectid, rat, d, e)
%function showConcurrent(selectid)
%function showConcurrent(selectid, rat, d, e)
% show concurrent firing across cell within a given epoch
% d :   day
% e :   epoch

% find tw-wide windows in which cells are co-firing
tw= 0.020; % [sec]
%tw= 2; % [sec]
twp= 0.5;  % [sec] display window width

setRoot;
if nargin<2
    epfile= sprintf('%s/data/epochs-%s', root, selectid);
    load(epfile);
    ep= epochs;
    nep= length(ep.rat);
else
    ep.rat{1}= rat;
    ep.num(1,:)= [d e];
    nep= 1;
end

cl= collectCellList(selectid);

for ie=1:nep
    rat= ep.rat{ie}; d= ep.num(ie,1); e= ep.num(ie,2);
    fprintf(1, '%s [%d %d]\n', rat,d,e);

    % load behavior data
    if ie==1 | ~strcmp(rat, ep.rat{ie-1}) | d~= ep.num(ie-1,1)
        bfile= sprintf('%s/%s/data2/behavdata%.2d.mat',root,rat,d);
        load(bfile);
        bd= behavdata{d}{e};
        t= bd.time;
        nt= length(t);
        load(sprintf('%s/%s/data2/spikedata%.2d.mat',root,rat,d));
        sd= spikedata{d}{e};
    end

    ind= find(strcmp(cl.rat, rat)' & cl.cellnum(:,1)==d & cl.cellnum(:,2)==e);

    ind= ind([6 4 1 2])  %@@@

    nc= length(ind);
    if nc<2; 
        warning('less than two cells in epoch'); 
        continue;
    end
    sp= {};
    for i=1:nc
        sp{i}= sd{cl.cellnum(ind(i),3)}{cl.cellnum(ind(i),4)};
        
%        spind= find(bd.traj(sp{i}.index)<0);
        spind= find(bd.ripple(sp{i}.index));
        sp{i}.time= sp{i}.time(spind);
        sp{i}.index= sp{i}.index(spind);
        sp{i}.linpos= sp{i}.linpos(spind);
        sp{i}.phase= sp{i}.phase(spind);
    end


    sw= round(tw/(t(2)-t(1)));
    swp= round(twp/(t(2)-t(1)));
    nw= ceil(nt/sw);
    count= zeros(1, nw);
    for ic=1:nc
        i= unique(floor(sp{ic}.index/sw)+1);
        count(i)= count(i)+1;
    end

    % # of co-firing  cells
    pair= find(count>=2);
    if isempty(pair); 
        warning('no co-occuring spikes found');
        continue;
    end

    % lump together neighboring theta cycles
    %hi= 1:length(pair)-1;  % show all windows 
    hi= find(diff(pair)>twp/tw);
    lo= sw*([ pair(1) pair(hi+1)]-1);
    hi= floor(sw*([pair(hi) pair(end)] +2));

    lo(lo<1)= 1;
    hi(hi>nt)= nt;
    fprintf(1, 'showing %d snippets\n', length(lo));

    % show all concurrent firing
    auxShowRasters(bd, sp, lo, hi);
end


function auxShowRasters(bd, sp, lo, hi)

figure;
set(gcf, 'Position', [6 549 1270 396]);
nc= length(sp);
nt= length(lo);

for it=1:nt
    ilo= lo(it); ihi= hi(it);
    tlo= bd.time(lo(it));
    thi= bd.time(hi(it));
    clf
    subplot(2,1,1)
    x= bd.linpos(ilo:ihi);

%    if any(x>20); continue; end %@@

    minx= min(x); maxx= max(x);
    if(maxx-minx<1); minx= minx-1; maxx=maxx+1; end
    plot(bd.time(ilo:ihi), x, 'k');  
    theta= bd.phase(ilo:ihi);
    theta= (1-cos(theta))/2*(maxx-minx)+minx;
    hold on
    plot(bd.time(ilo:ihi), theta, 'b');
%    axis([tlo thi 0 160]);
    set(gca, 'XLim', [tlo thi]);
    xlabel('time (s)');


    subplot(2,1,2)
    for ic=1:nc
        ind= find(tlo<sp{ic}.time & sp{ic}.time<thi);
        if isempty(ind); continue; end
        x= (sp{ic}.time(ind)-tlo)'*1000;
        nsp= length(x);
        lh= line([x; x], [ic-0.45; ic+0.45]*ones(1,nsp));
        set(lh, 'Color', 'r');
        set(lh, 'LineWidth', 1);
    end
    lh= line([0; (thi-tlo)*1000]*ones(1,floor(nc/2)), [1;1]*[2:2:nc]);
    set(lh, 'Color', 'k');
    set(lh, 'LineStyle', ':');
    axis ij
    axis([0 (thi-tlo)*1000 0.5 nc+0.5]);
    xlabel('\Delta t (ms)');
    str= input('(s)kip', 's');
    if(strcmp(str, 's')) break; end
end

