
create_filter = 0;
run_ff = 0;
save_ffdata = 0;
load_ffdata = 0;
stack_spikes = 0;
plotfigs = 0;

pilotEpoch = 1;
loadEpoch = 0;
make_licks = 0;
run_dfa = 1;

%% data filter params
Fp.animals = {'JZ1'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_licktrigspiking';
Fp.add_params = {'savefigs', 'wtrackdays', 'valid_ntrodes', 'exemplar_wepochs', ...
    'nonMU_cells'};
%% FF
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator,'eegtetrodes',Fp.tetfilter,...
        'cells', Fp.cellfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F); for d = 1:length(F); F(d).datafilter_params = Fp; end; end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s', Fp.epochEnvironment)); end
if load_ffdata
    spikesF = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s', Fp.epochEnvironment)); end

%% stack spikes
if stack_spikes
    spikestack = stack_riptrigspiking(spikesF);
end
%%
if plotfigs
    ian = 1;
    day = 4;
    epoch = 4;
    ntrode = 2;
    try
        suspikes = squeeze(spikestack(ian).suspikes{day}{epoch}(:,:,ntrode))';
    catch
        suspikes = [];
    end
    [xx, yy] = find(suspikes);
    s1 = scatterhist(xx/1000-1.001,yy, 'Kernel', 'on', 'bandwidth', ...
        [.02; 1], 'location','SouthEast', 'Direction', 'in', 'Marker', ...
        '+', 'color', 'k', 'MarkerSize', 2);
    % s1(3).Position(3) = .08;
    % s1(2).Position(4) = .08;
    axis tight
    ylabel('licknum','FontSize',8,'FontWeight','bold', 'FontName','Arial')
    xlim([-Fp.window(1) Fp.window(2)]); xticks([-Fp.window:.2:Fp.window]);
    xlabel('time s','FontSize',8,'FontWeight','bold', 'FontName','Arial');
    
    line([0 0],ylim, 'color','red', 'linewidth', 1)
    
end
%% pilot one epoch.. bypasses iterator
if pilotEpoch
    day = F(1).epochs{1}(1,1);
    epoch = F(1).epochs{1}(1,2);
    ntrode = F(1).data{1}{1}(end,1);
    cluster = F(1).data{1}{1}(end,2);
    animal = Fp.animals{1};
    animdef = animaldef(animal);
    if loadEpoch
        DIO = loaddatastruct(animdef{2}, animal, 'DIO', day);
        task = loaddatastruct(animdef{2}, animal, 'task', day);
        spikes = loaddatastruct(animdef{2}, animal, 'spikes', day);
        cellinfo = loaddatastruct(animdef{2}, animal, 'cellinfo', day);
    end
    %% make the licks events and save
    if make_licks
        get_licks(animal, [day epoch], DIO, task);
    else
        lick = load_data('filterframework', 'lick', animal, 'animpos', 0);
%         lick = loaddatastruct(animdef{2}, animal, 'lick');
    end
    if run_dfa
        epidx = find(ismember(F.epochs{1}, [day epoch], 'rows'));
        excludeIntervals = F.excludetime{1}{epidx};
        dfa_licktrigspiking([day epoch ntrode cluster], excludeIntervals, 'DIO', DIO, ...
            'task', task, 'spikes', spikes, 'cellinfo', cellinfo, 'events', ...
            lick, 'window', Fp.window, 'binsize', Fp.binsize);
    end
end



