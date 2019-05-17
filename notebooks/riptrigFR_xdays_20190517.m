

animals = {'JZ4'};
filtfunction = 'dfa_riptrigspiking';
filenamesave = 'riptrigFRxdays';
epenv = 'wtrack';

loaddata = 1;
plotfigs = 1;
pausefigs = 0;
savefigs = 1;

%% load data
if loaddata
    paths = make_paths(filtfunction, env);
    data = load_filter_output(paths.filtOutputDirectory, paths.filenamesave, ...
        animals);
    out = combine_riptrigspiking_perntrode(data);
end

%%
if plotfigs
    Pp = load_plotting_params(filenamesave);
    paths = make_paths(filenamesave, epenv);
    for ani = 1:length(animals)
        animal = animals{ani}
        ntrodes = [out{ani}.ntrode];
        for n = 1:length(ntrodes)
            %% ---- init fig----
            if savefigs && ~pausefigs;
                close all
                ifig = figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            %%
            nt_mu_keys = cell2mat({out{1}(n).mudata.index}');
            eplengths = cellfun(@(x) length(x(:,1)), {out{1}(n).mudata.instantFR}, 'un', 1);
            dayepvec = cell2mat(arrayfun(@(x,y,z) repmat([y z],x,1),eplengths',nt_mu_keys(:,1),nt_mu_keys(:,2),'un',0));
            a = cell2mat({out{1}(n).mudata.instantFR}');
            b = log(a+1); % adding one so that there's no negative log firing rates
            c = b(:,3:end-2);
            im = imagesc(c);
            colormap(jet)
            colorbar
            ylabel('swr#')
            xlabel('time')
            caxis([1 6])
            %% Super Axes
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %s %s nt%d', filenamesave, epenv, animal, ntrodes(n));
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', ...
                'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 14);
            
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(paths.figdirectory, paths.filenamesave, sprtit)
                close all
            end
        end
    end
end

if plotfigs
    Pp = load_plotting_params(filenamesave);
    paths = make_paths(filenamesave, epenv);
    for ani = 1:length(animals)
        animal = animals{ani}
        ntrodes = [out{ani}.ntrode];
        for n = 1:length(ntrodes)
            %% ---- init fig----
            if savefigs && ~pausefigs;
                close all
                ifig = figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            %%
            nt_mu_keys = cell2mat({out{1}(n).mudata.index}');
            eplengths = cellfun(@(x) length(x(:,1)), {out{1}(n).mudata.instantFR}, 'un', 1);
            dayepvec = cell2mat(arrayfun(@(x,y,z) repmat([y z],x,1),eplengths',nt_mu_keys(:,1),nt_mu_keys(:,2),'un',0));
            a = cell2mat({out{1}(n).mudata.instantFR}');
            b = log(a+1); % adding one so that there's no negative log firing rates
            c = b(:,3:end-2);
            im = imagesc(c);
            colormap(jet)
            colorbar
            ylabel('swr#')
            xlabel('time')
            caxis([1 6])
            %% Super Axes
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %s %s nt%d', filenamesave, epenv, animal, ntrodes(n));
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', ...
                'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 14);
            
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(paths.figdirectory, paths.filenamesave, sprtit)
                close all
            end
        end
    end
end