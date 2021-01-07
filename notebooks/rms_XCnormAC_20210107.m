

% load JZ4 day 1


% compute XCnormAC, measure RMS
% do permutation, measure RMS

% where is the observed RMS in the ranking of perm RMS distrubition

pconf = paramconfig;

create_filter = 0;
run_function = 1;
save_results = 0;
load_results = 0;

plot_perAn = 0;
plot_perTrial = 0;
show_figs = 0;
pause_figs = 0;
save_figs = 0;
save_fig_as = {'pdf', 'png'};


% specify data filters, iterators, functions into Fp struct
Fp.Label = 'rmsXCnormAC'; % analysis label (for results, plots)
Fp.animals = {'D10', 'D12', 'D13', 'JZ4'};
Fp.filtfunction = 'dfa_rmsXCnormAC';
Fp.params = {'firstToLastWellVisit', 'day1', 'wtrack', Fp.Label, ...
    Fp.filtfunction};
Fp = load_filter_params(Fp);

% run filter creation
if create_filter
    F = createfilter('animal', Fp.animals, 'days', Fp.days, 'epochs', ...
        Fp.epochfilter, 'excludetime', Fp.timefilter, 'iterator', ...
        Fp.iterator);
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