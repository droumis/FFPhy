%{

% plot position of ripples


%}



pconf = paramconfig;
create_filter = 1;
run_ff = 1;
load_ffdata = 0;

%% plot
plotfigs = 1;
showfigs = 1;
pausefigs = 0;
savefigs = 0;
savefigas = {'pdf', 'png'};

plot_ripPos_pAn_pDay = 1;
plot_ripPos_pAn = 0;
plot_ripPos = 0;
%% FF Data
Fp = [];
Fp.animals = {'JZ1'}; %{'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_ripPos'; % city.alien % not using space anymore
% expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
Fp.Label = 'ripPos';
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

%% Plot per day for each animal. then per animal. then all
if plotfigs
    if plot_ripPos_pAn_pDay
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
                
                plot(idata.txy(:,2), idata.txy(:,3), ...
                    'Color', [.8 .8 .8 .2], 'LineWidth', 1);
                hold on;
                axis tight
                xlimx = xlim;
                ylimy = ylim;
                colormap('jet')
                s=scatter(idata.SWR_xpos, idata.SWR_ypos, Pp.size, ...
                    idata.SWR_iLB, 'filled', 'MarkerFaceAlpha', .2);
%                 colorbar
                s.MarkerEdgeColor = 'none';
                xlim(xlimx)
                ylim(ylimy)
            end
            % super
            stit = sprintf('%s %s', figname, animal);
            setSuperAxTitle(stit);
            hold off
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' figname],...
                    stit, 'savefigas', savefigas);
            end
        end
    end

end