function showAF(selectid, adaptid, nsel, opt)
%function showAF(selectid, adaptid, nsel, opt)

maxt= 180; %@@
timeid= 'occ';
%maxt= 12;
%timeid= 'passes';

setRoot;
olddir= pwd;

load([root '/data/analist-' selectid])
nC= length(analist.rat);
if nargin<3; nsel=[1:nC]; end
tm= allTimes(timeid, selectid);
sp= allSpikeData(selectid, 1, nsel);
pf= allPlaceFields(selectid);

if nargin<4; 
    opt= [];
end

if ~isfield(opt, 'scale'); opt.scale= 1; end
if ~isfield(opt, 'putlabels'); opt.putlabels= 1; end
if ~isfield(opt, 'puttitle'); opt.puttitle= 1; end
if ~isfield(opt, 'units'); opt.units= 'rad'; end

switch opt.units
case 'rad'
    maxy= 2*pi;
    yfac= 1;
    ytick= [0 2 4 6];
case 'deg'
    maxy= 360;
    yfac= 180/pi;
    ytick= [0 180 360];
otherwise
    error('unknow unit');
end

Nsel= length(nsel);

oldrat= ''; oldnum= -1*ones(1,4);
npass= 0;
nspikes= 0;

global fmaux adaptest 

for iC= 1:nC
    isel= find(iC== nsel);
    if isempty(isel); continue; end

    rat= analist.rat{iC};
    num= analist.cellnum(iC,:); d=num(1); e=num(2); tet=num(3); c=num(4);
    if ~strcmp(oldrat, rat) 
        oldrat= rat; 
        cd([root '/' rat '/results/' adaptid]);
        setLocalOptions;
        fmaux.selectid= selectid;
        fmaux.select=[ fmaux.data2dir '/select-' fmaux.selectid];
        load(['stats-' selectid]);
        load(fmaux.select);
    end

    data=loadData(num);
    loadVar('.','adaptest',num(1));

    t= tm{iC};
    ind= find(isfinite(t));
    ind= [20 60 150 160]; %%@@
    opt.t= t(ind);
    opt.traj= analist.traj(iC);
%    opt.nx= 200; opt.ny= 70; % small print quality
    opt.nx= 600; opt.ny= 210; % large print quality

    out= eval2DModel(data,adaptest{d}{e}{tet}{c}.model, opt);
    x= out.x; y= out.y*yfac;


    iana= analist.iana(iC);
%    keyboard
%    statname= 'MutualInfo';
    statname= 'LinCorrOffset';
%    statstime= stats{d}{e}{tet}{c}.x{iana}.time;
    statstime= stats{d}{e}{tet}{c}.([statname '_time']){iana};
    svar= stats{d}{e}{tet}{c}.(statname){iana}(:,1);

%    selid= getSelectId([d e tet c]);
%    offset= select.a{selid}{1}.offset
    offset= unique(stats{d}{e}{tet}{c}.(statname){iana}(:,2))*yfac
    [tmp,cut]= min(mod(y+offset,maxy));
%    cut=1;

    [z{ind}]= deal(out.z{:}); 
    % find max
    if opt.scale
        zmax= 0;
        for iz=ind
            tmp= max(max(z{iz}));
            if tmp> zmax; zmax= tmp; end 
        end
        zmax
    end


    fh= figure;
    set(fh, 'Position', [357 54 604 892]);
    set(fh, 'Name', sprintf('%s [%d %d %d %d], no. %d', rat, d, e, tet, c, iC));
    nt= min(max(ind), maxt);
    padding= diff(pf(iC,:))/10;
    xmin= pf(iC,1)-padding;
    xmax= pf(iC,2)+padding;

%    for it=[55,155] %%@@
    for it=1:nt
        if isempty(z{it}); continue; end
%        subplot(ceil(nt/3), 3, it); %%@@
        figure

        ztmp= z{it}';
        ztmp= [ztmp(cut:end,:); ztmp(1:cut-1,:)];
        if opt.scale
            imagesc(x,y,ztmp, [0 zmax]);
            zmax
        else
            imagesc(x,y,ztmp);
            colorbar
        end
        axis xy
        if opt.putlabels
            xlabel('position (cm)');
            ylabel(['theta phase (' opt.units ')']);
        end

        set(gca, 'XTick', [100, 130]);
        set(gca, 'YTick', ytick);
%        set(gca, 'YTick', [0, pi]);
%        set(gca, 'YTickLabel', {'0' '3.14'});

        istats= min(find(t(it) < statstime));
%        if(abs(statstime(istats)-t(it)) > 5) error('could not find good time'); end

        if opt.puttitle
%            title(sprintf('pass %2d',it));
            ht= text(pf(iC,1), .5, sprintf('r= %.2f',abs(svar(istats))));
%            ht= text(pf(iC,1), .5, sprintf('%2d',it));
            set(ht, 'Color', 'w');
            set(ht, 'VerticalAlignment', 'bottom');
%            set(ht, 'HorizontalAlignment', 'right');
        end
%        title(sprintf('pass %d, t= %5.2f',it, t(it)));
%        title(sprintf('t= %5.2f, I=%.2f',t(it), I(istats)));
%        title(sprintf('%.2f',svar(istats)));
        fprintf(1, '%s %d, t= %5.2f, %s=%.2f\n', timeid, it, t(it), statname, svar(istats));
%        lh= line([1;1]*pf(iC,:), [0;4*pi]*[1,1]);
%        axis([xmin, xmax, 0, 4*pi]);
        axis([xmin, xmax, 0, maxy]);
%        lh= line([1;1]*pf(iC,:), [0;2*pi]*[1,1]);
%        set(lh, 'Color', 'w');
%        set(lh, 'LineWidth', 1);
%        axis square

        %@@
        figname= sprintf('%s/AF-%s-%.2d-%d%s', olddir, selectid, iC, it, timeid);
        myprint('mini', figname);
    end
%    orient tall
%    jointfig(gcf, ceil(nt/3), 3); %% @@
%    figname= sprintf('%s/AF-%s-%.2d', olddir, selectid, iC);
%    myprint(4.8*[1 4/3] , figname, 'paper-large', 0);

%    keyboard
end
cd(olddir);
fprintf(1, 'totals: %d spikes, %d passes, %.2f spikes/pass, median= %d\n', ...
    sum(nspikes), npass, mean(nspikes), median(nspikes));
