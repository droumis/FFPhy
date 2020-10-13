%{

XP-triggered average of z scored ripple power trace
as suggested by david
Meeting Notes: https://docs.google.com/document/d/1wTUtn20BHk7W64m5979wIZ19XYAZmgIi7sacT9YJoQc/edit#bookmark=id.fi6y0ntpmyqv

copy how i got z rip pwr from lickswrExamples_20190926
%}


pconf = paramconfig;
create_filter = 1;
run_ff = 1;
load_ffdata = 0;

%% plot
plotfigs = 0;
showfigs = 1;
pausefigs = 0;
savefigs = 1;
savefigas = {'pdf', 'png'};

plot_XPtrigAvgRip_pAn_pDay = 1;
plot_XPtrigAvgRip_pAn = 1;
plot_XPtrigAvgRip = 1;
%% FF Data
Fp = [];
Fp.animals = {'D10', 'D12', 'D13','JZ4'};
Fp.filtfunction = 'dfa_XPtrigAvgRip'; % city.alien % not using space anymore
% expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
Fp.Label = 'XPtrigAvgRip';
Fp.params = {'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell', ...
    'ripples>2', Fp.Label, Fp.filtfunction};

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
    if plot_XPtrigAvgRip_pAn_pDay
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
                % mean trace per an per day
                sf1 = subaxis(ndays, ncols, (d-1)*(ncols)+1, Pp.posparams{:});
                t = idata.time;
                m = idata.mean_rippwr_XPtrig;
                s = idata.sem_rippwr_XPtrig;
                
                st = knnsearch(t', Pp.win(1));
                en = knnsearch(t', Pp.win(2));
                t = t(st:en);
                m = m(st:en);
                s = s(st:en);
                
                plot(t, m, 'color', [0 0 0 1], 'linewidth', 1);
                hold on;
                fill([t'; flipud(t')],[m'-s';flipud(m'+s')], 'k',...
                    'linestyle', 'none', 'facealpha', .2);
                title(sprintf('day %d', d));
                ylabel('rippwr')
                xticks(Pp.win(1):.5:Pp.win(2))
                xlabel('time from XP');
                axis tight

            end
            % super
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
    if plot_XPtrigAvgRip_pAn
        figname = sprintf('%s-pAn',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        for a = 1:length(F)
            animal = F(a).animal{3};
            ifig = init_plot(showfigs, Pp.position);
            
            % plot            
            rippwr_XPtrig = vertcat(F(a).output{1}.rippwr_XPtrig);
            m = nanmean(rippwr_XPtrig,1);
            s = nanstd(rippwr_XPtrig,1)/sqrt(size(rippwr_XPtrig,1));
            t = idata.time;
            
            st = knnsearch(t', Pp.win(1));
            en = knnsearch(t', Pp.win(2));
            t = t(st:en);
            m = m(st:en);
            s = s(st:en);
            
            plot(t, m, 'color', [0 0 0 1], 'linewidth', 1);
            hold on;
            fill([t'; flipud(t')],[m'-s';flipud(m'+s')], 'k',...
                'linestyle', 'none', 'facealpha', .2);
            ylabel('rippwr')
            xticks(Pp.win(1):.5:Pp.win(2))
            xlabel('time from XP');
            axis tight;
            % super
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
    if plot_XPtrigAvgRip
        figname = sprintf('%s',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        
        % process
        rippwr_XPtrig = [];
        XPpostSWR_offset = [];
        for a = 1:length(F)
            rippwr_XPtrig = [rippwr_XPtrig; vertcat(F(a).output{1}.rippwr_XPtrig)];
        end
        m = nanmean(rippwr_XPtrig,1);
        s = nanstd(rippwr_XPtrig,1)/sqrt(size(rippwr_XPtrig,1));
        t = idata.time;
        st = knnsearch(t', Pp.win(1));
        en = knnsearch(t', Pp.win(2));
        t = t(st:en);
        m = m(st:en);
        s = s(st:en);
        % plot
        ifig = init_plot(showfigs, Pp.position);
        plot(t, m, 'color', [0 0 0 1], 'linewidth', 1);
        hold on;
        fill([t'; flipud(t')],[m'-s';flipud(m'+s')], 'k',...
            'linestyle', 'none', 'facealpha', .2);
        ylabel('rippwr')
        xticks(Pp.win(1):.5:Pp.win(2))
        xlabel('time from XP');
        axis tight
        % super
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