


pconf = paramconfig;
eventTrigSWR = 1; % NOW.. PIPE:space.city == per condition SWR mod (still need to somehow do per condition)

eventType = 'lick'; %lick swr
% run FF
create_filter = 0;
run_ff = 1;
load_ffdata = 0;

%% plot
plotfigs = 1;
showfigs = 1;
pausefigs = 0;
savefigs = 1;
savefigas = {'png', 'eps'};

plotSWRXPmod_perAn = 1;

%% FF Data
Fp = [];
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ3', 'JZ4'};
Fp.areas = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};

Fp.filtfunction = 'dfa_lickswrcorr'; % city.alien % not using space anymore
expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
Fp.Label = 'wXPTrigSWR';
Fp.params = {'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell', ...
    'ripples>2', Fp.Label, Fp.filtfunction};

Fp = load_filter_params(Fp);

if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'excludetime', ...
        Fp.timefilter, 'iterator', Fp.iterator);
    
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

%% Plot. per days for each animal. then per animal. then all
if plotfigs
    if plotSWRXPmod_perAn
        figname = 'wXPmodSWR';
        Pp=load_plotting_params({'defaults','dfa_lickswrcorr', figname});
        for a = 1:length(F)
            animal = F(a).animal{3};
            ifig = init_plot(showfigs, Pp.position);
            days = cell2mat([{F(a).output{1}.index}']);
            ndays = size(days,1);
            nc = 3;
            dNormxcZ = {};
            for d = 1:ndays
                day = days(d,1);
                idata = F(a).output{1}(d);
              %% Xcorr norm 
                nSWRs = length(idata.swrLickPhase);
                nSWRsN = 1-nSWRs/100;
                if nSWRsN < 0
                    nSWRsN = 0;
                end
                g = [nSWRsN nSWRsN nSWRsN];
                sf1 = subaxis(ndays+1, nc, (d-1)*(nc)+1, Pp.posparams{:});
                bar(idata.time, idata.normxc, 'FaceColor', g)
                xticks([])
%                 yticks([])
                ylim([0 .025])
                title(sprintf('%d Bouts %d XP. %d swr', size(idata.boutTimes, 1), length(idata.licks), nSWRs));
                ylabel(sprintf('%d', day));
              %% z-scored
                sf1 = subaxis(ndays+1, nc, (d-1)*(nc)+2, Pp.posparams{:});
                dNormxc = idata.smthxc;
                dNormxcSh = idata.smthxcSh;
                dNormxcZ{d} = (dNormxc - nanmean(dNormxcSh,1)) ./ nanstd(dNormxcSh,[],1);
                plot(idata.time, dNormxcZ{d}, 'color', g)
                xticks([])
%                 yticks([])
                ylim([-10 10])
                title(sprintf('%d Bouts %d XP. %d swr', size(idata.boutTimes, 1), length(idata.licks), length(idata.swrLickPhase)));
              %% phase
                sf2 = subaxis(ndays+1, nc, (d-1)*(nc)+3, Pp.posparams{:});
                polarhistogram(idata.swrLickPhase, 20, 'Normalization', 'pdf', 'facecolor', g)
%                 rticks([])
                thetaticks([])
                rlim([0 .8])
                pax = gca;
                pax.ThetaAxisUnits = 'radians';
                pax.ThetaZeroLocation = 'left';
                
            end
            %% all days normx
            subaxis(ndays+1, nc, (nc*(d-1))+4, Pp.posparams{:});
            swrXPxc = nanmean(cell2mat({F(a).output{1}.smthxc}'));
            bar(idata.time, swrXPxc, 'facecolor', 'b')
            ylabel('all');
            ylim([0 .02])
            %% all days z-scored
            NormxcZ = cat(1,dNormxcZ{:});
            NormxcZsem = sem(NormxcZ,1);
            NormxcZmean = mean(NormxcZ,1);
            
            subaxis(ndays+1, nc, (nc*(d-1))+5, Pp.posparams{:});
            plot(idata.time, NormxcZmean, 'b')
            hold on
            fill([idata.time'; flipud(idata.time')],[NormxcZmean'-NormxcZsem';flipud(NormxcZmean'+NormxcZsem')],'k',...
                'linestyle','none', 'facealpha', .2);
            ylim([-10 10])
             %% all days polar
            subaxis(ndays+1, nc, (nc*(d-1))+6, Pp.posparams{:});
            swrXPph = cell2mat({idata.swrLickPhase}');
            polarhistogram(swrXPph, 20, 'Normalization', 'pdf', 'facecolor', 'b')
            rlim([0 .8])
            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            pax.ThetaZeroLocation = 'left';
            %%
            stit = sprintf('%s %s daily', figname, animal);
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' figname], stit);
            end
        end
    end
end















