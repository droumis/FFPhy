
%{
What does JZ4 lfp pwr look like for the first X days
eventTrigLFP::
    - forest.bear.cactus.mushroom.beer.leaf
    // dfa_eventTrigLFP (bear:: ) ->
        -> stack_riptriglfp (cactus:: ) ->
        -> computeAnalyticSignal (mushroom:: ) ->
        -> makeEventSet (beer:: ) ->
        -> get_power (cactus:: eventSetPowerTrace) ->
- based on work from eventtrig_20191117.m
%}

pconf = paramconfig;
eventTrigLFP = 1; % PIPE:forest.bear.cactus.mushroom.beer.leaf == eventSet mean Spect
eventType = 'lick'; %lick swr

% run FF
create_filter = 0;
run_ff = 0;
load_ffdata = 1;

%% LFP to Analytic Signal Stack
% stack data
stack_LFP = 0; % Comely cactus
load_LFPstack = 0;
% get power, phase of all
make_rawpwr = 0; % Mendacious Mushroom
load_rawpwr = 0;

%% DESIGN MAT
make_expvarCat = 0; % Baleful beer
load_expvarCat = 0;

%% LFP Mean Power Per ExpVar
make_expvarCatMeanPwr = 0; % Lachrymose leaf
load_expvarCatMeanPwr = 0;
% combine per area
combineArea = 0;

%% LFP ITPC Per ExpVar

%% Spike Event Time Modulation
calcSUMod = 0; % Wheedling wheelbarrow
loadSUMod = 0;
gatherResults = 0; % combines across animals. keeps eventSet results seperate per unit

%% Spike Event Phasic Modulation


%% SWR x Event Phasic Modulation


%% plot
plotfigs = 0;
showfigs = 1;
pausefigs = 1;
savefigs = 1;

% swr
plotEventSU = 0;         % per eventSet, per SU
plotEventHeatRast = 0; % per eventSet, per area, per animal and all animals
plotEventHeatRastAllAni = 0; % requires gatherResults
plotEventModCDF = 0;     % requires gatherResults// per eventSet, per area, per animal and all animals

% LFP
plotLFPPerAreaAllAn = 1;

%%
Fp = [];
Fp.animals = {'D10'}; %, };
Fp.areas = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};

if eventTrigLFP
    Fp.filtfunction = 'dfa_eventTrigLFP'; % Bellicose Bear
    if strcmp(eventType, 'lick')
        expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
        Fp.Label = 'wtrackLickTrigLFP';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
        'excludeAfterLastWell', 'referenced', '4-350Hz',  Fp.Label, Fp.filtfunction};
    elseif strcmp(eventType, 'swr')
        expvars = {'all', 'lickbouts', 'nolickbouts'};
        Fp.Label = 'wtrackSWRTrigLFP';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
        'excludeAfterLastWell', 'referenced', '4-350Hz',  'ripples', ...
        Fp.Label, Fp.filtfunction};
    end
    
    
elseif eventTrigSpiking
    Fp.filtfunction = 'dfa_eventTrigSpiking'; % Redolent Rat
    if strcmp(eventType, 'lick')
        expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
        Fp.Label = 'wtrackLickTrigSpiking';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
            'excludeAfterLastWell', 'nonMU_cells', Fp.Label, Fp.filtfunction};
    elseif strcmp(eventType, 'swr') %'excludeNoise', 
        expvars = {'all', 'lickbouts', 'nolickbouts'};
        Fp.Label = 'wtrackSWRTrigSpiking';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
            'excludeAfterLastWell', 'nonMU_cells', 'ripples', ...
            Fp.Label, Fp.filtfunction}; % 'excludeNoise',
    end
end
Fp = load_filter_params(Fp);
if eventTrigLFP
    wp = getWaveParams(Fp.waveSet);
end
%% vectorize swr-trig lfp
if stack_LFP
    lfpstack = stack_riptriglfp(F, Fp);
end
if load_LFPstack
    lfpstack = load_data(Fp.paths.resultsDirectory, ['tensor_' Fp.Label], Fp.animals);
end

%% make rawpwr (all trials) [ntrode time rip freq]
if make_rawpwr
    % i need to fix this so that it doesn't crash anything but a super computer
    % make it memory aware, and calculate how many workers i can run at
    % once. 
    [rawpwr, ~] = computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'saveOutput',1, ...
        'lfptype', Fp.uselfptype, 'env', Fp.env, 'eventType', Fp.eventType); % uses parfor
end
if load_rawpwr
    rawpwr = load_data(sprintf('%s/analyticSignal/', pconf.andef{2}), ...
        sprintf('LFPpower_%s_%s_%s_%s', wp.waveSet, Fp.uselfptype, Fp.env, ...
        Fp.eventType), Fp.animals);
end
%% LFP POWER per condition
if make_expvarCatMeanPwr % :expvarCat @mean /ntTF $time
    evMPwr = getPower(expvarCat, rawpwr, Fp, 'run_perm', 0, 'eventType', Fp.eventType);
end
if load_expvarCatMeanPwr
    outdir = 'expvarCatMeanPwr';
    outpath = [pconf.andef{2},outdir,'/'];
    evMPwr = load_data(outpath, [outdir,'_', Fp.env '_' Fp.eventType], Fp.animals);
end
%% LFP combine perArea perCondition
if combineArea
    for ani = 1:length(evMPwr) % for each animal
        try
            animal = evMPwr(ani).animal{3};
        catch
            animal = evMPwr(ani).animal;
        end
        evMPwrArea(ani).animal = animal;
        evMPwrArea(ani).expvars = evMPwr(ani).expvars;
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
%         ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'') && 'ref'');
%         ntrodes = unique(ntrodes(:,3));
        ntrodes = evMPwr(ani).ntrode;
        for ia = 1:length(Fp.areas)
            ntsInArea = evaluatefilter(ntinfo,...
                sprintf('isequal($area,''%s'') && isequal($subarea,''%s'')', Fp.areas{ia}{1}, ...
                Fp.areas{ia}{2}));
            ntsInArea = unique(ntsInArea(:,3));
            ntsAIdx = knnsearch(ntrodes, ntsInArea);
            for iv = 1:length(evMPwr(ani).expvars)
                areaData = evMPwr(ani).data{iv}.pwr_mean_db(ntsAIdx,:,:);
                evMPwrArea(ani).data{ia,iv} = squeeze(nanmean(areaData,1))';
            end
        end
    end
%     evMPwrAreaAllAn = {evMPwrArea.data};
end
%% PLOT=====================================================================
if plotfigs
    
    %% plot LFP PWR per area per condition
    if plotLFPPerAreaAllAn
        if strcmp(eventType, 'swr')
            figname = 'wtrackSWRPwrAreaPerAn';
        elseif strcmp(eventType, 'lick')
            figname = 'wtrackLickPwrAreaPerAn';
        end
        Pp=load_plotting_params({'defaults',figname});
        for ani = 1:length(evMPwrArea) % for each animal
            animal = evMPwrArea(ani).animal;
            for ia = 1:length(Fp.areas)
                ifig = init_plot(showfigs, Pp.position);
                for iv = 1:length(evMPwrArea(ani).expvars)
                    sf = subaxis(1,length(Fp.areas), iv);
                    frequency = round(wp.frex);
                    time = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
                    trim = knnsearch(time',Pp.win(1)):knnsearch(time', Pp.win(2));
                    ttime = time(trim);
                    tdata = evMPwrArea(ani).data{ia,iv}(:,trim);
                    contourf(ttime,frequency,tdata,40,'linecolor','none')
                    set(gca,'ydir','normal','yscale','log');
                    
                    colormap(Pp.usecolormap)
                    if iv == 1
                        caxis(sf, 'auto')
                        cax = [-max(abs(caxis)) max(abs(caxis))];
                    end
                    caxis(sf, cax) % normalize to the 'all' condition
                    
                    
%                     if ia == 1
%                         caxis([-4 4]); %cax(2)])
%                     else
%                         caxis([-1.5 1.5]); %cax(2)])
%                     end
                    ytickskip = 2:4:wp.numfrex;
                    set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', Pp.tickFsize)
                    title(sprintf('%s',evMPwr(ani).expvars{iv}))
                    yl = ylim;
                    xl = xlim;
                    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                    colorbar
                    hold on
                end
%                 ylabel(sf{1},'frequency Hz')
%                 xlabel(sf{1}, 'time s')
                %%
                stit = sprintf('%s %s %s %s %s %s %s',figname, Fp.areas{ia}{1},Fp.areas{ia}{2}, Fp.eventType, ...
                    animal, Fp.env, Fp.waveSet);
                setSuperAxTitle(stit);
                if pausefigs
                    pause
                end
                if savefigs
                    save_figure([pconf.andef{4} figname], stit);
                end
            end
            
        end
    end
end