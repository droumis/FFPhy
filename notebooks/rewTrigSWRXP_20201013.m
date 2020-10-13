%{


reward triggered XP and SWR 

%}

pconf = paramconfig;
create_filter = 0;
run_ff = 0;
load_ffdata = 0;

%% plot
plotfigs = 1;
showfigs = 1;
pausefigs = 0;
savefigs = 0;
savefigas = {'pdf', 'png'};

plot_pAn_pDay = 0;
plot_pAn = 1;
plot_all = 0;
%% FF Data
Fp = [];
Fp.animals = {'JZ1'}; %{'D10', 'D12', 'D13','JZ4'};
Fp.filtfunction = 'dfa_rewTrigSWRXP'; % city.alien % not using space anymore
% expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
Fp.Label = 'rewTrigSWRXP';
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

%% Plot per day for each animal. then per animal. then all
if plotfigs
    if plot_pAn_pDay
        figname = sprintf('%s-pAn-pDay',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        for a = 1:length(F)
            animal = F(a).animal{3};
            ifig = init_plot(showfigs, Pp.position);
            days = cell2mat([{F(a).output{1}.index}']);
            ndays = size(days,1);
            ncols = 1;
            for d = 1:ndays
                day = days(d,1);
                idata = F(a).output{1}(d);
                sf1 = subaxis(ndays, ncols, (d-1)*(ncols)+1, Pp.posparams{:});
                %% PLOT
                
            end
            %% super
            stit = sprintf('%s %s', figname, animal);
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' figname],...
                    stit, 'savefigas', savefigas);
            end
        end
    end
    if plot_pAn
        figname = sprintf('%s-pAn',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        for a = 1:length(F)
            animal = F(a).animal{3};
            ifig = init_plot(showfigs, Pp.position);
            %% swr raster
            idata = F(a).output{1};
            time = idata(1).time';
            sT = knnsearch(time, -abs(Pp.win(1)));
            eT = knnsearch(time, Pp.win(2));
            
            sf = subaxis(6, 1, [1 2], Pp.posparams{:});
            rewTrigSWR = cell2mat({idata.rewTrigSWR}');
            [xx,yy] = find(rewTrigSWR(:,sT:eT)');
            swrTwin = (xx*1e-3)-abs(Pp.win(1));
            h = scatter(swrTwin, yy, Pp.spikeSz, '.k', 'markeredgealpha', ...
                Pp.spikeAlpha);
            axis tight
            xlim([-abs(Pp.win(1)) Pp.win(2)])
%             xticks([]);
            ylabel('SWR #');
            line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5], 'linewidth', 2)
            
            %% swr PSTH
            sf = subaxis(6,1,3,Pp.posparams{:});
            [xx,yy] = find(rewTrigSWR');
            swrTwin = (xx*1e-3)-abs(time(1));
            bintime = [time(1):Pp.bin:time(end)]';
            sB = knnsearch(bintime, -abs(Pp.win(1)));
            eB = knnsearch(bintime, Pp.win(2));
            h = histc(swrTwin, bintime);
%             hs = smoothdata(h, 1,'loess', 4);
            area(bintime(sB:eB), h(sB:eB), 'facecolor', 'k')
            axis tight
            xlabel('time s')
            ylabel('count')
            line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5], 'linewidth', 2)
            
            %% xp raster
            idata = F(a).output{1};
            time = idata(1).time';
            sT = knnsearch(time, -abs(Pp.win(1)));
            eT = knnsearch(time, Pp.win(2));
            
            sf = subaxis(6, 1, [4 5], Pp.posparams{:});
            rewTrigXP = cell2mat({idata.rewTrigXP}');
            [xx,yy] = find(rewTrigXP(:,sT:eT)');
            swrTwin = (xx*1e-3)-abs(Pp.win(1));
            h = scatter(swrTwin, yy, Pp.spikeSz, '.k', 'markeredgealpha', ...
                Pp.spikeAlpha);
            axis tight
            xlim([-abs(Pp.win(1)) Pp.win(2)])
%             xticks([]);
            ylabel('SWR #');
            line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5], 'linewidth', 2)
            
            %% xp PSTH
            sf = subaxis(6,1,6,Pp.posparams{:});
            [xx,yy] = find(rewTrigXP');
            swrTwin = (xx*1e-3)-abs(time(1));
            bintime = [time(1):Pp.bin:time(end)]';
            sB = knnsearch(bintime, -abs(Pp.win(1)));
            eB = knnsearch(bintime, Pp.win(2));
            h = histc(swrTwin, bintime);
%             hs = smoothdata(h, 1,'loess', 4);
            area(bintime(sB:eB), h(sB:eB), 'facecolor', 'k')
            axis tight
            xlabel('time s')
            ylabel('count')
            line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5], 'linewidth', 2)
            
            %%
            allAxesInFigure = findall(gcf,'type','axes');
            linkaxes(allAxesInFigure, 'x');

            %% super
            stit = sprintf('%s %s', figname, animal);
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' figname],...
                    stit, 'savefigas', savefigas);
            end
        end
    end
    if plot_all
        figname = sprintf('%s',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        %% PLOT
        
        %% super
        stit = sprintf('%s', figname);
        setSuperAxTitle(stit);
        if pausefigs
            pause
        end
        if savefigs
            strsave = save_figure([pconf.andef{4} '/' figname],...
                stit, 'savefigas', savefigas);
        end
    end
end