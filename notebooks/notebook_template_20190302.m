

Fp = load_filter_params({'wtrack','ripples','nonref_ntrodes',...
    '>100spikes_cells','riptrigspiking'});
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ3', 'JZ4'};
% Fp.days = [1:7];

runFilterFramework = 1;
loadFilterOutput = 0;
plotstuff = 0;

%% ---------------- Paths ---------------------------------------------------
paths = make_paths(Fp);
%% ---------------- Run FIlter -----------------------------------------------
if runFilterFramework == 1    
    F = createfilter('animal',Fp.animals,'epochs', ... 
        Fp.epochfilter, 'cells',Fp.cellfilter, 'excludetime', Fp.timefilter, ...
        'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, {'spikes', ...
        Fp.eventDataLabel ,'pos','task'}, 'TF',Fp.TF,'window',Fp.window, ...
        'binsize',Fp.binsize,'frbinsize',Fp.frbinsize,'minthresh', ...
        Fp.minstdthresh,'maxvelocity',Fp.maxvelocity,'minvelocity', ...
        Fp.minvelocity, 'consensus_numtets',Fp.consensus_numtets,'welldist', ...
        Fp.welldist);
    F = runfilter(F);
    for a = 1:length(F) % save filter detailes along with results
        F(a).datafilter_params = Fp;
    end
end
%% ---------------- Save Filter Output ----------------------------------------
if saveFilterOutput == 1;
    save_filter_output(F, paths.filtOutputDirectory, paths.filenamesave)
end
%% ---------------- Load Filter Output ----------------------------------------
if loadFilterOutput == 1;
    F = load_filter_output(paths.filtOutputDirectory, paths.filenamesave);
end
%% ---------------- Plot ------------------------------------------------------
if plotstuff
    Pp = load_plotting_params('riptrigspiking');
    %% set loop, etc

    
    %% init new figure 
    if savefigs && ~pausefigs;
        ifig = figure('Visible','off','units','normalized','position', ...
            Pp.position);
    else
        ifig = figure('units','normalized','position',Pp.position);
    end
    set(gcf,'color','white')
    %% plot script
    
    
    

    
    %% ---- super title and colorbar----
    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    sprtit = sprintf('%s - %d %d %d', animalID, idata.dtc(1),...
        idata.dtc(2), idata.dtc(3));
    iStitle = text(.5, .95, {sprtit}, 'Parent', sprtitleax, 'Units', ...
        'normalized');
    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
        'horizontalAlignment', 'center');
     %% ---- pause, save figs ----
    if pausefigs
        pause
    end
    if savefigs
        save_figure(paths.figdirectory, paths.filenamesave, sprtit)
    end
    close all
end

