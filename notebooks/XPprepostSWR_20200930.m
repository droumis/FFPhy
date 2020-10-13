

%{

David wants histograms of SWR offsets from nearest XP, both pre and post
Meeting Notes: https://docs.google.com/document/d/1wTUtn20BHk7W64m5979wIZ19XYAZmgIi7sacT9YJoQc/edit#bookmark=id.fi6y0ntpmyqv

%}


% for each animal/day, load filtered rips, licks and find offset
% xpswrAll_20191105.m
% swrlickxcorr_20191031.m

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

plot_XPprepostSWR_pAn_pDay = 1;
plot_XPprepostSWR_pAn = 1;
plot_XPprepostSWR = 1;
%% FF Data
Fp = [];
Fp.animals = {'D10', 'D12', 'D13', 'JZ1','JZ4'};
Fp.filtfunction = 'dfa_XPprepostSWR'; % city.alien % not using space anymore
% expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
Fp.Label = 'XPprepostSWR';
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
    if plot_XPprepostSWR_pAn_pDay
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
                %% hist
                sf1 = subaxis(ndays, ncols, (d-1)*(ncols)+1, Pp.posparams{:});
                histogram(idata.XPpreSWR_offset_shuf, Pp.hbins, ...
                    'FaceColor', [.7 .7 .7], 'EdgeColor', [1 1 1]);
                hold on;
                histogram(idata.XPpostSWR_offset_shuf, Pp.hbins, ...
                    'FaceColor', [.7 .7 .7], 'EdgeColor', [1 1 1]);
                
                histogram(idata.XPpreSWR_offset, Pp.hbins, ...
                    'FaceColor', Pp.preClr)
                hold on;
                histogram(idata.XPpostSWR_offset, Pp.hbins, ...
                    'FaceColor', Pp.postClr)
                title(sprintf('day%d', d));
                xlabel('time from SWR');
                xlim([-.5 .5])
            end
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
    if plot_XPprepostSWR_pAn
        figname = sprintf('%s-pAn',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        for a = 1:length(F)
            animal = F(a).animal{3};
            ifig = init_plot(showfigs, Pp.position);
            %% hist
            %             sf1 = subaxis(1, ncols, (d-1)*(ncols)+1, Pp.posparams{:});
            XPpreSWR_offset_shuf = vertcat(F(a).output{1}.XPpreSWR_offset_shuf);
            XPpostSWR_offset_shuf = vertcat(F(a).output{1}.XPpostSWR_offset_shuf);
            histogram(XPpreSWR_offset_shuf, Pp.hbins, 'FaceColor', [.7 .7 .7], ...
            'EdgeColor', [1 1 1]);
            hold on;
            histogram(XPpostSWR_offset_shuf, Pp.hbins, 'FaceColor', [.7 .7 .7], ...
            'EdgeColor', [1 1 1]);
        
            XPpreSWR_offset = vertcat(F(a).output{1}.XPpreSWR_offset);
            histogram(XPpreSWR_offset, Pp.hbins, ...
                'FaceColor',Pp.preClr)
            hold on;
            XPpostSWR_offset = vertcat(F(a).output{1}.XPpostSWR_offset);
            histogram(XPpostSWR_offset, Pp.hbins, ...
                'FaceColor',Pp.postClr)
            xlabel('time from SWR');
            ylabel('XP count');
            xlim([-.5 .5])
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
    if plot_XPprepostSWR
        figname = sprintf('%s',Fp.Label);
        Pp=load_plotting_params({'defaults',Fp.Label, figname});
        XPpreSWR_offset = [];
        XPpostSWR_offset = [];
        XPpreSWR_offset_shuf = [];
        XPpostSWR_offset_shuf = [];
        for a = 1:length(F)
            XPpreSWR_offset = [XPpreSWR_offset; vertcat(F(a).output{1}.XPpreSWR_offset)];
            XPpostSWR_offset = [XPpostSWR_offset; vertcat(F(a).output{1}.XPpostSWR_offset)];
            
            XPpreSWR_offset_shuf = [XPpreSWR_offset_shuf; vertcat(F(a).output{1}.XPpreSWR_offset_shuf)];
            XPpostSWR_offset_shuf = [XPpostSWR_offset_shuf; vertcat(F(a).output{1}.XPpostSWR_offset_shuf)];
        end
        %% hist
        %             sf1 = subaxis(1, ncols, (d-1)*(ncols)+1, Pp.posparams{:});
        ifig = init_plot(showfigs, Pp.position);
        histogram(XPpreSWR_offset_shuf, Pp.hbins, 'FaceColor', [.7 .7 .7], ...
            'EdgeColor', [1 1 1]);
        hold on;
        histogram(XPpostSWR_offset_shuf, Pp.hbins, 'FaceColor', [.7 .7 .7], ...
            'EdgeColor', [1 1 1]);
        
        histogram(XPpreSWR_offset, Pp.hbins, 'FaceColor',Pp.preClr)
        hold on;
        histogram(XPpostSWR_offset, Pp.hbins, 'FaceColor',Pp.postClr)
        xlabel('time from SWR');
        ylabel('XP count');
        xlim([-.5 .5])
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