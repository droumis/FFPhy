
%{

what pct of ripples are during lick bouts?

- make a ven diagram of iLB, eLB 
- make a spatial plot of where ripples occur

%}


pconf = paramconfig;
create_filter = 1;
run_ff = 1;
load_ffdata = 0;

%% plot
plotfigs = 1;
showfigs = 1;
pausefigs = 0;
savefigs = 1;
savefigas = {'pdf', 'png'};

plot_pctILB_pAn_pDay = 0;
plot_pctILB_pAn = 0;
plot_pctILB = 1;

%% FF
Fp = [];
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ4'};
Fp.Label = 'pctILB';
Fp.filtfunction = 'dfa_pctILB';
Fp.params = {'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell', ...
    'ripples>2', 'proximalWell', Fp.Label, Fp.filtfunction};

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

%% plot

if plotfigs
    if plot_pctILB_pAn_pDay
        figname = sprintf('%s-pAn-pDay',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        ifig = init_plot(showfigs, Pp.position);
        for a = 1:length(F)
            animals{a} = F(a).animal{3};
            days = cell2mat([{F(a).output{1}.index}']);
%             ndays = size(days,1);
            ncols = 1;
%             for d = 1:ndays
%             day = days(d,1);
%             idata = F(a).output{1}(d);
            
            sf1 = subaxis(1, ncols, (d-1)*(ncols)+1, Pp.posparams{:});
            numSWRiLB = vertcat(F(a).output{1}.numSWRiLB);
            numSWR = vertcat(F(a).output{1}.numSWR);
            pctSWRiLB = numSWRiLB ./ numSWR * 100;
            plot(1:length(pctSWRiLB), pctSWRiLB, 'linewidth', 3)
            hold on
%                 title(sprintf('day%d', d));
%             end
        end
        legend(animals)
        xlabel('day');
        ylabel('pct iLB');
        title('>2 std')
        %%
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
    if plot_pctILB_pAn
        figname = sprintf('%s-pAn',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        for a = 1:length(F)
            animal = F(a).animal{3};
            ifig = init_plot(showfigs, Pp.position);
            %% hist
            %             sf1 = subaxis(1, ncols, (d-1)*(ncols)+1, Pp.posparams{:});
            XPpreSWR_offset = vertcat(F(a).output{1}.XPpreSWR_offset);
            histogram(XPpreSWR_offset, 20)
            hold on;
            XPpostSWR_offset = vertcat(F(a).output{1}.XPpostSWR_offset);
            histogram(XPpostSWR_offset, 20)
            xlabel('time from SWR');
            ylabel('XP count');
            xlim([-.25 .25])
            %%
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
    if plot_pctILB
        figname = sprintf('%s',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        numSWRiLB = [];
        numSWR = [];
        for a = 1:length(F)
            numSWRiLB(a,:) = nansum(vertcat(F(a).output{1}.numSWRiLB));
            numSWR(a,:) = nansum(vertcat(F(a).output{1}.numSWR));
            animals{a} = F(a).animal{3};
        end
        pctSWRiLB = numSWRiLB ./ numSWR * 100;
        %% hist
        %             sf1 = subaxis(1, ncols, (d-1)*(ncols)+1, Pp.posparams{:});
        ifig = init_plot(showfigs, Pp.position);
        plot(repmat([1:size(pctSWRiLB,2)], size(pctSWRiLB,1), 1)', pctSWRiLB', 'linewidth', 2)
        hold on
        d = reshape(repmat([1:size(pctSWRiLB,2)], size(pctSWRiLB,1), 1), ...
            [numel(pctSWRiLB),1]);
        scatter(d, pctSWRiLB(:), 1000, '.k')
%         notBoxPlot(pctSWRiLB)
        b = boxplot(pctSWRiLB, 'Symbol', '.','Color', 'k');
        set(gca,'XTickLabel',{' '})
        xticks([1 2 3 4])
        xticklabels({'2', '3', '4', '5'})
        
%         set(gca, 'FontSize', 10)
%         legend off
%         ylim([0 50]);
        xlabel('ripple std threshold');
        ylabel('pct iLB');
        legend(animals, 'Location','northwest')
        %%
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