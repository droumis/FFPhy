%{


reward triggered XP and SWR 

%}

pconf = paramconfig;
create_filter = 1;
run_ff = 1;
load_ffdata = 0;

%% plot
plotfigs = 0;
showfigs = 0;
pausefigs = 0;
savefigs = 0;
savefigas = {'pdf', 'png'};

plot_pAn_pDay = 0;
plot_pAn = 0;
plot_all = 0;
%% FF Data
Fp = [];
Fp.animals = {'JZ1'}; %{'D10', 'D12', 'D13','JZ4'};
Fp.filtfunction = 'dfa_rewTrigSWRXP'; % city.alien % not using space anymore
% expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
Fp.Label = 'rewTrigSWRXP';
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
            %% PLOT


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