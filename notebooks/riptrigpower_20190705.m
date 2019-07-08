


run_ff = 0;
savedata = run_ff;


load_stack = 0;
get_ripstate = 0;
calc_AS = 0;
batch_load_as = 0;
combine_and_subsample = 0;
calc_PWR = 0;
load_PWR = 1;
run_permtest = 0;

plotfigs = 0;
savefigs = 0;
pausefigs = 0;

% use_filters = {'firstwell', 'noise'};
Fp.animals = {'D10'};%, 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
Fp.add_params = {'wtrack', 'wavelets4-300Hz', 'excludeNoise','excludePriorFirstWell'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
me = animaldef('Demetris');

filetail = '';
%% (ntrode x time x ripple)
if load_stack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end

%% Calculate and Save Analytic Signal
if calc_AS
    computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'overwrite', 1);
end

%% get behavioral state for each rip
if get_ripstate
    ripstate = getStateFilters(lfpstack);
end

%% load analytic signal
% i should really seperate out the conditions into different files
% parfor'd this takes 6 min all 30 or 12 seconds per ntrode to load 1 lfptype from an animal. 80 Gb
% i don't this works with the parfor in getpower.. maybe bc it's trying to
% allocate the full as loaded to each worker?
% without parfor, each ntrode takes 60 seconds to load so 30 minutes in
% loading time

if batch_load_as
    tic
    lfptypes = {'eeg'};
    asAnim = cell(1,length(Fp.animals));
    for ani = 1:length(Fp.animals)
        animal = Fp.animals{ani};
        anidx = find(strcmp(animal, {lfpstack.animal}));
        ntrodes = lfpstack(anidx).ntrodes;
        as = cell(1,length(lfptypes));
        for itype = 1:length(lfptypes)
            asnt = cell(1,length(ntrodes));
            parfor nti = 1:length(ntrodes)
                nt = ntrodes(nti);
                asnt{nti} = loadAS(Fp.animals{anidx}, nt, Fp.waveSet, lfptypes{itype});
            end
            as{itype} = asnt;
            clear asnt
        end
        asAnim{ani} = as;
        clear as
    end
    batchtime = toc
end
%% combine all ntrode data from D10, condition 'all' into a single matrix
if combine_and_subsample
    tic
    a = [asAnim{1}{1}{:}];
    b = {a.as};
    asdatacat = cat(4,b{:});
    % subsample the data
%     asdatacatsub = asdatacat(1:2:end,:,:,:);
    cmb_time = toc
end
%% Compute and Save Power from analytic signal for each ntrode, condition
if calc_PWR
    % without parfor, .4 minutes for power + 10 permutations, 1 ntrode
    % (.4* 30) * (100*30) = 36000 minutes = 600 hours or 25 days for 1000 permutations
    % now .2 per so 25/2 = 12.5 days without parfor
    % 10 permutations takes 15 seconds.. so 1000 would take 1500 seconds or
    % 25 minutes per ntrode.. 25*30 = 750 minutes or 12.5 hours total for all 1000
    % permutes for all ntrodes 12.5 hours + (80-15) * 30/60/60 = 13 hours
    % as things are now.. if i can find a way to parfor the ntrodes.. what
    % if instead i combined them in a matrix? ok i matrixafied it.. 
    % ok so now it's down to 11.5 minutes with 10 perms all nts
    % so 1000 perms would be about 19 hours. maybe a lot less
    % 100 perms took about 90 minuts so 1000 would actually probably take
    % 10.5 hours.. which is an overnight.. not bad.. subsampling to half
    % the time would still take several hours, which is still an
    % overnight.. unless there's another 10x.. might as well keep it the
    % way it is.. but to do 7 animals in a reasonable time.. i do need
    % another 10x.. especially if i want to run more conditions, lfptypes,
    % measures.. 
    % typhoon has a shitty 750 GEforce gpu so there's no way i can fit on
    % it's free mem with gpuArray
    % - now trying to downsample by 2.. no parfor (was breaking even with
    % as few as 8 workers).. i could try even just 2 workers.. 
    % 50 seconds for compute_power: 2 permutes, 15 ntrodes ,downsampled by 2
    % 50*500*2/60/60
    % - 33 sec: 1 perm, all ntrodes, downs 100
    % 1.7 s: 1perm, 1 nt, dsamp 10 1.7*1000*30/60/60 = 14 hours
    % 18 s: 100perm, 1nt, dsamp10 18*10*30/60/60 = 1.5 hours
    % 180s: 1000perm, 1nt, dsamp10  180*30/60/60 = 1.5 hours
    % 104s: above, with parfor2  104*30/60 = 52 min
    % 46s: above, with parfor4  46*30/60 = 23 min
    % 54s: parfor2 1perm 30nt dsamp10  54*1000/60/60 = 15 hours
    % 52s: same as above but partfor6  
    % 83s: parf4 10perm 30nt dsamp10  83*100/60/60 = 2.3 hours  * this was weird.. zmap took 70 s 
    % 310s : parf4 10perm 30nt dsamp2   8 hours * dsamped earlier at the AS
    % use :parf4 1000perm 30nt dsamp10 .. should take 2.3 hours for D10
    %   all wtrack power + permtest
    % run then plot all ntrodes mean power + zmap overlay
    % then run the D12,D13,JZ1. i still need to do JZ2-4 asignal i think
    % then plot all the animals ntrodes
%     tic
    getPower(asdatacat, lfpstack, ripstate, Fp, 'lfptypes', {'eeg'}, ...
        'ripstatetypes', {'all'});
%     pwr_time = toc
end

%% plot per ntrode

%% combine by area, across animals, ntrodes
%% plot per area (ca1, supdors, supvent, deepdors, deepvent)
%% 
%{
to validate the ntrode locations across ntrodes.. i need to look at the
power/phase results of each ntrode

i also need to validate the ca1 position before combining ca1 tets

sup first 

first figure 
    - recording, wtrackbehavior
    - riptrig power spect of each of area
    - mean lfp trace for each spect
    - what else?

second figure is like first but showing ITPC

third figure - PWRCORR between all areas

fourth figure - ISPC between all areas

fifth figure - power, powercorr, itpc, ispc for correct, incorrect, diff

sixth figure - power, powercorr, itpc, ispc for outbound, inbound, diff

seventh figure - power, powercorr, itpc, ispc for awake, sleep, diff

eighth figure - power, powercorr, itpc, ispc for performance, perf change

nineth figure - power, powercorr, itpc, ispc for epoch phase, day phase



%}





