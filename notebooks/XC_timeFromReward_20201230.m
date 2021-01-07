
%{
get xcnormac as a function of time since reward

need: reward output time via dio
need: 


%}

create_filter = 1;
run_function = 1;
save_results = 1;
load_results = 0;

plot_perAn = 0;
plot_perTrial = 0;
show_figs = 1;
pause_figs = 0;
save_figs = 0;
save_fig_as = {'pdf', 'png'};

% specify data filters, iterators, functions into Fp struct
Fp.Label = 'XCfromRew'; % analysis label (for results, plots)
Fp.animals = {'JZ4'};
Fp.filtfunction = 'dfa_XCfromRew';
Fp.params = {'firstToLastWellVisit', 'day1', 'wtrack', Fp.Label, Fp.filtfunction};
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
        
        d = vertcat(F(a).output{1}.xcNormAc_incr);
        for incr = 1:size(d,2)
%         semd = std(d,[],2)/sqrt(size(d,2));
%         t = -.5:1/1500:.5;
%         dm = nanmean(d,2);
%         fill([t'; flipud(t')],[dm-semd;flipud(dm+semd)], 'b',...
%             'linestyle', 'none', 'facealpha', .2);
%         hold on
            plot(t, dm, 'color', rand(1,3), 'linewidth', 2, 'DisplayName', ...
                sprintf('%d from rew', incr-1));
            ylim([-.5 .5])
            line([0 0], ylim, 'linestyle', '--', 'color', 'k')
            xlabel('time from XP')
        end
        legend
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







