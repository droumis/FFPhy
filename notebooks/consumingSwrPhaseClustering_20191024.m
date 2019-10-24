
%{
yesterday i collected ILI phase for swrs...
today i need to
--- create column of consumption lickbout vs nonconsumption
i think an initial first pass would be to label the swr with whether it
occurred within ~5 seconds of reward output.. what does the distribution of
lick burst duration look like?
- for each animal
    - load the DIO times, task
- for each swr
    - grab the Dout ID from task and collect those Douts' times from DIO
- lookup swr start time from nearest preceding reward out
- if yes, add 1 to column for that swr

- then do a kuiper test (circ ks test) to compare distributions
%}
maxTimeSinceRew = 5;
pconf = paramconfig;

create_filter = 0;
run_ff = 0;
save_ffdata = 0;
load_ffdata = 0;
createResultStruct = 0;
getConsumeCol = 0;
calcWetVDry = 0;

plotfigs = 1;
savefigs = 1;
pausefigs = 0;

%% get the lickphase for each swr
% data filter params
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_lickswrcorr';
Fp.params = {'savefigs', 'wtrackdays', Fp.filtfunction};
% FF
Fp = load_filter_params(Fp);
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F);
    for d = 1:length(F); F(d).datafilter_params = Fp; end
end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', sprintf('_%s', Fp.epochEnvironment));
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s', Fp.epochEnvironment));
end

%% collect ILI phase results
if createResultStruct
    for ani = 1:length(F)
        if 0
            andef = animaldef(F(ani).animal);
            ripInfo = loaddatastruct(andef{2}, F(ani).animal, 'ca1rippleskons');
        end
        Fs = cell2mat(F(ani).output');
        swrLickPhase = cell2mat({Fs.swrLickPhase}');
        swrStart = cell2mat({Fs.swrInBurstStart}');
        swrBurstInterval = cell2mat({Fs.swrBurstInterval}');
        deUnq = cell2mat({Fs.index}');
        de = [];
        for i = 1:size(deUnq)
            de = [de; repmat(deUnq(i,:), length(Fs(i).swrInBurstStart), 1)];
        end
        Os(ani).data = struct;
        Os(ani).animal = F(ani).animal;
        for s = 1:length(swrStart)
            riptimes = ripInfo{de(s,1)}{de(s,2)}{1}.starttime;
            sIdx = knnsearch(riptimes, swrStart(s));
            Os(ani).data(s).day = de(s,1);
            Os(ani).data(s).epoch = de(s,2);
            Os(ani).data(s).swrLickPhase = swrLickPhase(s);
            Os(ani).data(s).starttime = ripInfo{de(s,1)}{de(s,2)}{1}.starttime(sIdx);
            Os(ani).data(s).endtime = ripInfo{de(s,1)}{de(s,2)}{1}.endtime(sIdx);
            Os(ani).data(s).duration = Os(ani).data(s).endtime - Os(ani).data(s).starttime;
            Os(ani).data(s).maxthresh = ripInfo{de(s,1)}{de(s,2)}{1}.maxthresh(sIdx);
            Os(ani).data(s).swrBurstStart = swrBurstInterval(s,1);
            Os(ani).data(s).swrBurstEnd = swrBurstInterval(s,2);
            %     Os(s).starttimecheck0 = swrStart(s) - Os(s).starttime;
        end
    end
end
%% for each swr, get whether it was near consumption or not (wet v dry)
if getConsumeCol
    for ani = 1:length(Os)
        andef = animaldef(Os(ani).animal);
        task = loaddatastruct(andef{2}, Os(ani).animal, 'task');
        DIO = loaddatastruct(andef{2}, Os(ani).animal, 'DIO');
        animal = Os(ani).animal;
        
        for x = 1:length(Os(ani).data)
            day = Os(ani).data(x).day;
            epoch = Os(ani).data(x).epoch;
            
            % from lick DIO ID in taskInfo to index into DIO
            rewDIOID = task{day}{epoch}.outputdio;
            isOutput = cellfun(@(x) isequal(x.input,0), DIO{day}{epoch}, 'un', 1);
            dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')),...
                DIO{day}{epoch}, 'un', 1);
            rewDIOIdx = find(all([ismember(dioID,rewDIOID)' isOutput'],2));
            rewPortDOut = DIO{day}{epoch}(rewDIOIdx);
            
            % for each lick port Din, gather it's times along with index
            l = cellfun(@(x) x.times, rewPortDOut,'un',0)';
            [rewTime,X] = sort(cell2mat(l));
            g = [];
            for d = 1:length(rewDIOID)
                g{d,1} = repmat(rewDIOID(d), length(l{d}),1);
            end
            i = cell2mat(g);
            id = i(X);
            
            e = Os(ani).data(x).starttime - rewTime;
            [timeSinceRew, idxNearestRew] = min(e(e > 0));
            
            Os(ani).data(x).timeSinceRew = timeSinceRew;
            Os(ani).data(x).wet = timeSinceRew < maxTimeSinceRew;
        end
    end
end

%% Plot the two distributions
if plotfigs
for ani = 1:length(Os)
    animal = Os(ani).animal;
    wVdSWR = [Os(ani).data.swrLickPhase; Os(ani).data.wet]';
    figname = 'wetVdryILIphaseSWR';
    wetSWR = wVdSWR(wVdSWR(:,2)==1,1);
    drySWR = wVdSWR(wVdSWR(:,2)==0,1);
    
    Pp=load_plotting_params({'defaults',figname}, 'savefigs', savefigs);
    polarhistogram(wetSWR, 32, 'Normalization', 'pdf', 'edgealpha', .1, 'FaceColor', 'b');
    hold on
    polarhistogram(drySWR, 32, 'Normalization', 'pdf', 'edgealpha', .1, 'FaceColor', 'r');
%     polarhistogram([Os(ani).data.swrLickPhase]', 32, 'Normalization', 'pdf', 'edgealpha', .1, 'FaceColor', 'k');
    legend({'wet', 'dry'})
    hold off
    [pval, k, K] = circ_kuipertest(wetSWR, drySWR, 100, 0);
    sprintf('p:%d', pval)
    title([figname sprintf(' p:%d', pval)])
    if savefigs
        save_figure('/stelmo/demetris/analysis/figures/wetVdry/', [figname animal])
    end
end
end

%% 
% why is this different than what i had before for swr illi phases?? i've
% changed some things to invclude more swrs.. and i recoded the getting of
% the phase code.. i need to check this... 














