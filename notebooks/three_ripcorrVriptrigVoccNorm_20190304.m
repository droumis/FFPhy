
%{
- fix results saving to make it one per animal
- remake ripcorr, riptrig, behave, ratemap and save per animal
- make jz2, d12 riptrig spiking
- remake riptrig spiking and include all epochs
- compute corrcoef per swr in new analysis function
- load ripcorr, riptrig, behave, ratemap

- combine results per day, subplots tight, perform with conf bounds
- next: compute per rip measure
- next: compute autocorr w score
%}

animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
me = animaldef('demetris');
env = 'wtrack';
loadresults = 0;
combine_results = 1;
save_combines_results = 0;
compute_across_results = 0;
plotfigs = 0;

if load_results
    for an = 1:length(animals);
        f = 'mua_calcxcorrmeasures';
        ripcorr{an} = load(sprintf('%s/%s/%s_%s_%s.mat',me{2}, f, f, env, animals{an}));
        f = 'dfa_riptrigspiking';
        riptrig{an} = load(sprintf('%s/%s/%s_%s_%s.mat',me{2}, f, f, env, animals{an}));
        f = 'dfa_occNormFiring';
        ratemaps{an} = load(sprintf('%s/%s/%s_%s_%s.mat',me{2}, f, f, env, animals{an}));
        f = 'behave_state';
        performan{an} = load(sprintf('%s/%s/%s_%s_%s.mat',me{2}, f, f, env, animals{an}));
    end
end

if combine_results
    % combine the results for each day
    % should i make a new combined result or 1 for each result type ?
end

if compute_across_results
    % compute corr coef between rip corr and performance/change
end

% if plotfigs
%     %% animal loop
%     % for each animal, create a graph for the corr of the rip corr measure Vs ...
%     % performance, per day
%     % could also compute other across-days scores for rip trig, ratemaps,
%     % autocorr, etc
%     
%     %% figure (nt-pair) loop
% end
