% swr's are necessary for learning
% ripple chains in ca1 are more common on novel, and under memory load

% Detect, plot rip chains

% are these ripple chains that are more common early w track actually doing the back and
% forth with MEC?

% do more ntrodes, cells, mu fr rate, ripple power, fast gamm amplitude
% become active in latter the chainned?
% Measure MEC from chain to chain.. does it build up after the first?

% tonegawa: After detecting candidate ripple events (see Candidate Ripple and Replay Event Detection), time lags between peaks of ripple event were computed. Singlet ripples were defined as ripple events that are temporally separated by 200 ms or more to adjacent ripple events. For the doublets and triplets (ripple bursts), the temporal lags were set to less than 200 ms but greater than 70 ms. These parameters were determined by computing the single ripple duration (∼80 ms) and inter-ripple intervals (∼125 ms) when ripple bursts occurred. In the event that adjacent ripple closer than 70 ms, these were categorized as single ripple. Number of adjacent ripple peaks was determined and assigned to doublets, triplets and others (more peaks than three). We performed this analysis to both quiet awake and slow-wave sleep periods that are identified by elevated delta wave power (1-4 Hz), lowered theta wave power

run_ff = 0;
savedata = run_ff;
load_lfp = 0;
stack_lfp = 0;
save_spikes = 1;
run_spikes_ff = 1;

%
load_events = 0;
make_filter_vecs = 0;
get_chains = 0;
filterpureWdays = 0;

plotfigs = 0;
% plot_heatrast_traces_perep_allnt = 1;

savefigs = 1;
pausefigs = 0;

animals = {'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};

add_params = {'wtrack'};
Fp.animals = animals;
Fp.add_params = add_params;
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
filetail = '';
use_filters = {'firstwell', 'noise'};
%% run filter/func
if run_ff
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    F = runfilter(F);
    for d = 1:length(F)
        F(d).datafilter_params = Fp;
    end
end
%% save data
if savedata == 1
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s%s', Fp.epochEnvironment, filetail))
end
%% load lfp
if load_lfp
    F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s%s', Fp.epochEnvironment, filetail));
end
%% stack LFP
if stack_lfp
    lfpstack = stack_riptriglfp(F);
end
%%
if 0
    for ani = 1:length(lfpstack)
        days = unique(lfpstack(ani).day);
        numchained = {};
        for d = 1:length(days)
            daystarts = lfpstack(ani).ripStartTime(lfpstack(ani).day == day);
            dayends = lfpstack(ani).ripEndTime(lfpstack(ani).day == day);
            day = days(d);
            t = 1;
            numr = length(daystarts);
            chain = zeros(numr,1);
            for i = 2:numr
                interv = daystarts(i) - dayends(i-1);
                if interv > .02 && interv < .5
                    chain(i-1) = t;
                    chain(i) = t;
                else
                    chain(i) = 0;
                    t=t+1;
                end
            end
            
            lfpstack(ani).chaingroup{d,1} = chain;
            % get [numsinglets, numdoublets, ...]
            lfpstack(ani).numchained{d}(1) = sum(chain == 0);
            unqchrip = unique(chain);
            unqchrip = unqchrip(2:end); % get rid of the zero
            for i = 1:length(unqrip)
                chainlen(i) = sum(chain == unqrip(i));
            end
            for j = 2:max(chainlen)
                lfpstack(ani).numchained{d}(j) = sum(chainlen==j);
            end
        end
    end
end
%%
% since the lfop is based on ca1 swrs.. i can't use it for detecting mec
% rips.. but i can just load the mec rips from the ripkons and speed, std
% threshold it (could even use the exclude times way)..

% show how ca1 chains change over days
% show how mec chains change over days
% very few.. need to lower std threshold.. load the ripples from ripkons
% instead of the riptrig lfp.. 

% get activation of mec at the each  how mec chains change over days
% what would be interesting either way it comes out is to see whether the
% same ntrodes/cells respond in chains, or if the response changes from one
% rip to the next in terms of the strength of response at different ntrode
% sites.. 

% get the ripple chains from th ripkons

%% run spiking FF
sFp.animals = animals;
add_params = {'wtrack'};
sFp.add_params = add_params;
sFp.filtfunction = 'dfa_riptrigspiking';
sFp = load_filter_params(sFp, 'add_params',sFp.add_params);

if run_spikes_ff == 1
    spikesF = createfilter('animal',sFp.animals, 'epochs', sFp.epochfilter, ...
        'cells',sFp.cellfilter, 'excludetime',sFp.timefilter, 'iterator',sFp.iterator);
    spikesF = setfilterfunction(spikesF,sFp.filtfunction,sFp.datatypes,sFp.options{:});
    spikesF = runfilter(spikesF);
    for d = 1:length(spikesF)
        spikesF(d).datafilter_params = sFp;
    end
end
%% save spikes data
if save_spikes == 1
    save_data(spikesF,sFp.paths.filtOutputDirectory,sFp.paths.filenamesave, 'filetail',...
        ['_',sFp.epochEnvironment])
end
%% load spikes
if load_spikes
    spikesF = load_filter_output(sFp.paths.filtOutputDirectory,sFp.paths.filenamesave, ...
        sFp.animals, 'filetail', ['_',sFp.epochEnvironment]);
end

%% stack spikes
% if stack_spikes
%     spikestack = stack_riptrigspiking(spikesF);
% end


%% load
if get_chains
    for ani = 1:length(animals)
        animal = animals{ani};
        andef = animaldef(animal);
        
        % load events
        events(ani).animal = animal;
        events(ani).ca1ripples = loaddatastruct(andef{2}, animal, 'ca1rippleskons');
        events(ani).mecripples = loaddatastruct(andef{2}, animal, 'mecrippleskons');
        
        % filter events by velocity, std
%         pos(ani).pos = loaddatastruct(andef{2}, animal, 'pos');
        
        % load task
        task(ani).animal = animal;
        task(ani).info = loaddatastruct(andef{2}, animal, 'task');
        dayeps = evaluatefilter(task(ani).info, 'strfind($environment, ''track'')');
        if filterpureWdays
            % save list of ONLY pure wtrack days and epochs (not rotated, 6arm)
            envs = cellfetch(task(ani).info, 'environment');
            notwtrack = evaluatefilter(task(ani).info, ...
            '(~isequal($environment, ''openfield'') && ~isequal($environment, ''wtrack'') && isequal($type, ''run''))');
            dayeps = dayeps(~ismember(dayeps(:,1), notwtrack(:,1)),:);
        else
            dayeps = dayeps;
        end
        
        % get ripples filtered for std, speed, etc
        events(ani).ca1_rips = getconstimes(andef{2}, animal, dayeps,'ca1rippleskons', 1, ...
            'consensus_numtets',1, 'minstdthresh', 2,'exclusion_dur',0, ...
                'minvelocity', 0,'maxvelocity',4);
        
        events(ani).mec_rips = getconstimes(andef{2}, animal, dayeps,'mecrippleskons', 1, ...
            'consensus_numtets',1, 'minstdthresh', 2,'exclusion_dur',0, ...
            'minvelocity', 0,'maxvelocity',4);
        
        % get numspikes per mec su during (and maybe later within 50ms following
        % each ca1 swr)
        % label each swr whether it is the first second or third in a row
        % (includes singlets)
        
        % label each swr with which mec ntrode had the most # of spikes
        % 1 *** group the numspikes and plot the result as the magnitude response
        % 2 *** for the not singlets, does the max numspikes mec ntrode
        % match the preceding rip in that chain?
        
        for de = 1:length(dayeps(:,1))
            day = dayeps(de,1);
            epoch = dayeps(de,2);
            starttimes = events(ani).ca1_rips{day}{epoch}.starttime;
            endtimes = events(ani).ca1_rips{day}{epoch}.endtime;
            
            
        end
        
    end
end


























