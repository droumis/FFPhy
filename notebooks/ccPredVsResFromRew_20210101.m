%{

determine the corrcoef of XCnormAc-predicted ripple power trace vs the
real, as a function of time since reward 

%}
pconf = paramconfig;

create_filter = 0;
run_function = 0;
save_results = 1;
load_results = 0;

plot_perAn = 1;
plot_perTrial = 0;
show_figs = 1;
pause_figs = 0;
save_figs = 1;
save_fig_as = {'pdf', 'png'};


% specify data filters, iterators, functions into Fp struct
Fp.Label = 'ccPvRXCfromRew'; % analysis label (for results, plots)
Fp.animals = {'D10', 'D12', 'D13', 'JZ4'};
Fp.filtfunction = 'dfa_ccPvRXCfromRew';
Fp.params = {'firstToLastWellVisit', 'wtrackdays', Fp.Label, Fp.filtfunction};
Fp = load_filter_params(Fp);

% run filter creation
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, ...
        'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end

% run function
if run_function
   F = arrayfun(@(x) setfield(F(x), 'datafilter_params', Fp), 1:length(F), ...
       'un', 1); % save the datafilter params along with the animal data
   F = runfilter(F);
   if save_results
       save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
           'filetail', ['_', Fp.Label]);
   end
end

% load function results -> 
if load_results
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', ['_' Fp.Label]);
end

if plot_perAn
    figname = sprintf('%s-pAn',Fp.Label);
    Pp=load_plotting_params({'defaults',Fp.Label, figname});
    for a = 1:length(F)
        animal = F(a).animal{3};
        ifig = init_plot(show_figs, Pp.position);
        
        % get data
        coefs = vertcat(F(a).output{1}.coefs);
        t = 0:size(coefs, 2)-1;
        % plot
        dm = nanmean(coefs)';
        semd = nanstd(coefs)/sqrt(size(coefs,1))';
        fill([t'; flipud(t')],[dm-semd;flipud(dm+semd)], [.5 .5 .5],...
            'linestyle', 'none');
        hold on
        plot(t, dm, 'color', 'k', 'linewidth', 2);
        xlabel('time from rew (sec)')
        ylabel('corrcoef')
        
        % super axis
        stit = sprintf('%s %s', figname, animal);
        setSuperAxTitle(stit);
        if pause_figs
            pause
        end
        if save_figs
            strsave = save_figure([pconf.andef{4} '/' figname],...
                stit, 'savefigas', save_fig_as);
        end
    end
end
if plot_perTrial
    figname = sprintf('%s-pTrial',Fp.Label);
    Pp=load_plotting_params({'defaults',Fp.Label, figname});
    for a = 1:length(F)
        animal = F(a).animal{3};
        ifig = init_plot(show_figs, Pp.position);
        days = cell2mat([{F(a).output{1}.index}']);
        % get data per day
        for d = 1:size(days,1)
            day = days(d,1);
            idata = F(a).output{1}(d);
            predRipStack = idata.predRipStack;
            realRipStack = idata.realRipStack;
            rewIntervs = idata.rewIntervs;
            xp = idata.xp
            for iT = 1:size(rewIntervs,1)
                t = rewIntervs(iT,1):1500:rewIntervs(iT,2);
                % plot
                plot(t, predRipStack(iT,:), 'k', 'linewidth', 2);
                hold on;
                plot(t, realRipStack(iT,:), 'b', 'linewidth', 2);
                axis tight
                xpTr = isIncluded(xp, rewIntervs(iT,:));
                line([xpTr xpTr], ylim, 'color', [.5 0 0 .5]);
                xlabel('time from rew (sec)')
                
                % super axis
                stit = sprintf('%s %s %.03f', figname, animal, ...
                    rewIntervs(iT,1));
                setSuperAxTitle(stit);
                if pause_figs
                    pause
                end
                if save_figs
                    strsave = save_figure([pconf.andef{4} '/' figname],...
                        stit, 'savefigas', save_fig_as);
                end
            end
        end
    end
end






















