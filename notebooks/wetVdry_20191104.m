
%{
refactor of consumingSwrPhaseClustering_20191024.m
wet (rewarded/consuming) v dry (error) Burst-Lick phase clustering and rxn diff
DR 20191104
%}

pconf = paramconfig;
% get data
create_filter = 0;
run_ff = 0;
load_ffdata = 0;
% proc results
createResultStruct = 0;
getConsumeCol = 0;
% plot
plotfigs = 1;
pausefigs = 0;
savefigs = 1;
%%
Fp.animals = {'D10', 'D13', 'JZ1', 'JZ4'}; %, 'JZ4'}; %;{'D10'}; %
Fp.filtfunction = 'dfa_lickswrcorr';
Fp.Label = 'wetVdryLickSWRITPC';
Fp.params = {'wtrackdays', 'ripples', 'lickbouts', Fp.filtfunction};
Fp = load_filter_params(Fp);

%% FF
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.Label]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.Label]);
end
%%
if createResultStruct
    D = struct;
    for i = 1:length(F)
        Fs = cell2mat(F(i).output');
        D(i).animal = Fs(i).animal;
        D(i).swrLickPhase = cell2mat({Fs.swrLickPhase}');
        D(i).swrStart = cell2mat({Fs.swrInBurstStart}');
        D(i).swrBurstInterval = cell2mat({Fs.swrBurstInterval}');
        D(i).dayEpoch = cell2mat({Fs.dayEpoch}');
        D(i).dayEpochUniq = cell2mat({Fs.index}');
    end
end

%% for each swr, get whether it was near consumption or not (wet v dry)
% need to change this so that it's proximity of the start of the containing
% lick-burst rather than the time of each swr
if getConsumeCol
    for a = 1:length(D)
        animal = D(a).animal;
        andef = animaldef(animal);
        task = loaddatastruct(andef{2}, animal, 'task');
        DIO = loaddatastruct(andef{2}, animal, 'DIO');
        for x = 1:length(D(a).swrStart) % for each swr
            day = D(a).dayEpoch(x,1);
            ep = D(a).dayEpoch(x,2);
            % get reward DIO ID from taskInfo to index into DIO
            rewDIOID = task{day}{ep}.outputdio;
            isOutput = cellfun(@(x) isequal(x.input,0), DIO{day}{ep}, 'un', 1);
            dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')),...
                DIO{day}{ep}, 'un', 1);
            rewDIOIdx = find(all([ismember(dioID,rewDIOID)' isOutput'],2));
            rewPortDOut = DIO{day}{ep}(rewDIOIdx);
            
            % for each lick port Din, gather it's times along with index
            l = cellfun(@(x) x.times, rewPortDOut,'un',0)';
            [rewTime,X] = sort(cell2mat(l));
            g = [];
            for d = 1:length(rewDIOID)
                g{d,1} = repmat(rewDIOID(d), length(l{d}),1);
            end
            i = cell2mat(g);
            id = i(X);
            e = D(a).swrStart(x) - rewTime;
            timeSinceRew = min(e(e > 0));
            
            D(a).timeSinceRew(x,1) = timeSinceRew;
            D(a).isWetLick(x,1) = timeSinceRew < Fp.maxTimeSinceRew;
        end
    end
end
%% Plot the two distributions
if plotfigs
    figname = 'wetVdryILIphaseSWR';
    for a = 1:length(D)
        Pp=load_plotting_params({'defaults',figname}, 'savefigs', savefigs, 'pausefigs', pausefigs);
        subaxis(1,1,1,Pp.posparams{:})
        animal = D(a).animal;
        wVdSWR = [D(a).swrLickPhase D(a).isWetLick];
        wetSWR = wVdSWR(wVdSWR(:,2)==1,1);
        drySWR = wVdSWR(wVdSWR(:,2)==0,1);
        dh = polarhistogram(drySWR, 24, 'Normalization', 'pdf', 'edgealpha', .1, 'FaceColor', 'r');
        hold on
        wh = polarhistogram(wetSWR, 24, 'Normalization', 'pdf', 'edgealpha', .1, 'FaceColor', 'b');
        legend([wh dh], {'wet', 'dry'})
        hold off
%         [pval, k, K] = circ_kuipertest(wetSWR, drySWR, 100, 1);
        %%
        stit = sprintf('%s %s %s', figname, animal, Fp.env);     
        title(stit);
%         setSuperAxTitle([stit sprintf(' p:%.03f', pval)]);
        if pausefigs
            pause
        end
        if savefigs
            save_figure([pconf.andef{4} figname], stit);
        end
    end
end

                % smooth polar histogram??
%         n = 25;
%         % start fit
%         val = h.Values;
%         c = fftshift(fft(ifftshift(val)))/n;   % fourier coefficients
%         % n0 is index for constant term.  c(n0) = npts/n = average bin value
%         n0 = (n+1)/2;    
%         B1 = 2*abs(c(n0+2));
%         B2 = angle(c(n0+2));
%         theta1 = (h.BinEdges(1:end-1) + h.BinEdges(2:end))/2;   % bin centers
%         E = c(n0) + B1.*cos(2*theta1 + B2);
%         hold on
%         polarplot(theta1,E,'-o')
%         hold off
