

%{

XPtrigSWR

%}
create_filter = 0;
run_ff = 0;
load_ffdata = 0;

plotfigs = 1;
showfigs = 0;
pausefigs = 0;
savefigs = 1;
savefigas = {'pdf', 'png'};

plot_XPtrigSWR_pAn = 1;
% plot_XPtrigSWR = 0;
plot_XPtrigSWR_1SWR_pAn = 0;
%% FF Data
Fp = [];
Fp.animals = {'D10', 'D12', 'D13', 'JZ1'}; %'JZ4'
Fp.filtfunction = 'dfa_XPtrigSWR'; % city.alien % not using space anymore
% expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
Fp.Label = 'XPtrigSWR';
Fp.params = {'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell', ...
    'ripples>2', 'wetLickBursts', 'proximalWell', Fp.Label, Fp.filtfunction};

Fp = load_filter_params(Fp);

if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,...
        'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F),...
        'un', 1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.Label]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', ['_' Fp.Label]);
end



%% PLOT
if plotfigs
    if plot_XPtrigSWR_1SWR_pAn
        figname = sprintf('%s-pAn','XPtrigSWR-1SWR');
        Pp=load_plotting_params({'defaults',figname});
        for a = 1:length(Fp.animals)
            animal = Fp.animals{a};
            ifig = init_plot(showfigs, Pp.position);
            idata = F(a).output{1};
            time = idata(2).time';
            SWRtimeFromNearestXP_Sh = cell2mat({idata.SWRtimeFromNearestXP_Sh}');
            SWRtimeFromNearestXP_Sh = SWRtimeFromNearestXP_Sh(SWRtimeFromNearestXP_Sh<.1);
            SWRtimeFromNearestXP_Sh = SWRtimeFromNearestXP_Sh(SWRtimeFromNearestXP_Sh>-.1);
            histogram(SWRtimeFromNearestXP_Sh, Pp.hbins, ...
                'FaceColor', [.5 .5 .5], 'EdgeColor', [.5 .5 .5]);
            hold on;
            
            SWRtimeFromNearestXP = cell2mat({idata.SWRtimeFromNearestXP}');
            SWRtimeFromNearestXP = SWRtimeFromNearestXP(SWRtimeFromNearestXP>-.1);
            SWRtimeFromNearestXP = SWRtimeFromNearestXP(SWRtimeFromNearestXP<.1);
            histogram(SWRtimeFromNearestXP, Pp.hbins, ...
                'FaceColor', [.4 .1 .4], 'EdgeColor', [.4 .1 .4]);
            line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5],...
                'linewidth', 2)
            xlim([-abs(Pp.win(1)) Pp.win(2)])
            ylabel('SWR count');
            xlabel('SWR time from nearest XP')
                %% super
            stit = sprintf('%s %s', figname, animal);
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure('demetris', stit, 'savefigas', ...
                    savefigas, 'subdir', figname);
            end
        end
    end
    if plot_XPtrigSWR_pAn
        figname = sprintf('%s-pAn',Fp.Label);
        Pp=load_plotting_params({'defaults',figname});
        for a = 1:length(Fp.animals)
            animal = Fp.animals{a};
            ifig = init_plot(showfigs, Pp.position);
%             sf1 = subaxis(1, 3, 1, Pp.posparams{:}, 'SpacingHoriz', Pp.SpHz);
            %% swr raster
            subaxis(4, 1, [1 2], Pp.posparams{:});
            idata = F(a).output{1};
            time = idata(2).time';
            XPTrigSWR = cell2mat({idata.XPiLBtrigSWR_raster}');
            [xx,yy] = find(XPTrigSWR');
            xpTwin = (xx*1e-3)-abs(time(1));
            scatter(xpTwin, yy, Pp.ptSz, '.k', 'markeredgealpha', ...
                Pp.ptAlpha);
            axis tight
            xlim([-abs(Pp.win(1)) Pp.win(2)])
%             xticks([]);
            ylabel('XP #');
            line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5],...
                'linewidth', 2)
%             title('XPTrigSWR')
            
            %% swr PSTH
            subaxis(4,1,3,Pp.posparams{:});
            histogram(xpTwin, Pp.hbins, ...
                'FaceColor', [0 0 0], 'EdgeColor', [.5 .5 .5]);
            xlim([-abs(Pp.win(1)) Pp.win(2)])
            ylabel('count');
            %% swr avg rate
            subaxis(4,1,4,Pp.posparams{:});
%             [xx,yy] = find(XPTrigSWR');
%             swrTwin = (xx*1e-3)-abs(time(1));
            bintime = [time(1):Pp.bin:time(end)]';

            [hcell,hedges] = arrayfun(@(r) ...
                histcounts((find(XPTrigSWR(r,:))*1e-3)-abs(time(1)), bintime), ...
                (1:size(XPTrigSWR,1)).', 'Uni',0);
            hmtx = cell2mat(hcell); % Recover Numeric Matrix From Cell Array
            hmtx = hmtx * 1/Pp.bin; % put into Hz
            edges = hedges{1};
            centers = edges(2:end)-(diff(edges(1:2))/2);
            m = nanmean(hmtx)';
            hm = smoothdata(m, 1,'loess', 10);
            sem = (nanstd(hmtx)/sqrt(size(hmtx,1)))';
            hsem = smoothdata(sem, 1,'loess', 10);
            fill([centers; flipud(centers)],[hm-hsem;flipud(hm+hsem)], ...
                [.5 .5 .5], 'linestyle', ...
                'none', 'facealpha', .5);
            hold on;
            plot(centers, hm, 'color', [0 0 0 1], 'linewidth', 1)
            
            axis tight
            xlim([-abs(Pp.win(1)) Pp.win(2)])
            
            
%             [counts, centers] = hist(rewTrigXP');
%             hs = smoothdata(h, 1,'loess', 4);
%             area(bintime(sB:eB), h(sB:eB), 'facecolor', 'k')
%             [f,xi] = ksdensity(swrTwin, );
%             plot(xi,f);
%             xlim([-abs(Pp.win(1)) Pp.win(2)])
%             xlabel('time s')
            ylabel('Hz')
            line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5], 'linewidth', 2)

            
            %%
%             allAxesInFigure = findall(gcf,'type','axes');
%             linkaxes(allAxesInFigure, 'x');

            %% super
            stit = sprintf('%s %s', figname, animal);
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure('demetris', stit, 'savefigas', ...
                    savefigas, 'subdir', figname);
            end
        end
    end
end