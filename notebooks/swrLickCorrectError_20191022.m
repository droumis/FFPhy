

%{
*** update.. for each swr phase, also get the task novelty, performance, and learning rate for
each swr, alongside correct/error and duration.. so 1 matrix of Nswrs x 6
features..
*** maybe also do time from epoch/day start? yes
*** then do stats..
- for correct v error.. i need to compare the non-uniformity (Raleigh Z value?) of the two
distributions.
- for correct v error.. also maybe just do a circular version of the ks
test to just ask if the distributions are different from each other?,, this
would account for different phase preferences... even if there's no
difference in the non-uniformity.
- for novely, performance, learning-rate.. these are all continuous
covariates, so i need to do a circular-linear correlation

swr locking to lick phase (continuous) during correct vs error trials

get measure per animal

---- get lick phase of each swr

- lickbout intervals
- swr times (within lickbout)
- lick times (within lickbout)
=== all swr lick phases

---- group correct vs error
group into whether they were during error trials or correct trials
determine this by whether there was a rew output within ~2 seconds of
lickbout start
- rew DOuput times
=== correct, error swr lick phases

---- test for significance
shuffle the group identities 1k times to get a ctrl distribution.
could also do a ttest (testing that the means came from difference
distributions)

+++++++++ q's ++++++++++++
- are there enough error trial swrs to do this?
- how to combine animals?
- also get a measure of lickbout duration and report the correct vs error
duration distributions..
- does the result change if the correct lickbouts are randomly truncated to
match the distribution of the error bouts?
- is there a difference in lick phase locking within a lick bout? head v tail.
%}

%% get the lickphase for each swr
pconf = paramconfig;
create_filter = 0;
run_ff = 1;
save_ffdata = 1;
load_ffdata = 0;
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
if 0
    andef = animaldef(Fp.animals{ani});
    ripInfo = loaddatastruct(andef{2}, Fp.animals{ani}, 'ca1rippleskons');
end

%%
if 1
    %% Get index, Lick Phase of swr
    ani = 1;
    Fs = cell2mat(F(ani).output');
    swrLickPhase = cell2mat({Fs.swrLickPhase}');
    swrStart = cell2mat({Fs.swrInBurstStart}');
    swrBurstInterval = cell2mat({Fs.swrBurstInterval}');
    deUnq = cell2mat({Fs.index}');
    de = [];
    for i = 1:size(deUnq)
        de = [de; repmat(deUnq(i,:), length(Fs(i).swrInBurstStart), 1)];
    end
    
    %% Test phase modulation using all swrs
    [pval, z] = circ_rtest(swrLickPhase); % pval is the stat rayleigh test. z is mean res vec
    phasemod = log(z); % i think log makes it 'variance normalized' (karalis,sirota)
end
if 1
    %% lookup swr in ripplekons based on [day epoch starttime] dayepIdx swrStart
    Os = struct;
    for s = 1:length(swrStart)
        riptimes = ripInfo{de(s,1)}{de(s,2)}{1}.starttime;
        sIdx = knnsearch(riptimes, swrStart(s));
        Os(s).day = de(s,1);
        Os(s).epoch = de(s,2);
        Os(s).swrLickPhase = swrLickPhase(s);
        Os(s).starttime = ripInfo{de(s,1)}{de(s,2)}{1}.starttime(sIdx);
        Os(s).endtime = ripInfo{de(s,1)}{de(s,2)}{1}.endtime(sIdx);
        Os(s).duration = Os(s).endtime - Os(s).starttime;
        Os(s).maxthresh = ripInfo{de(s,1)}{de(s,2)}{1}.maxthresh(sIdx);
        Os(s).swrBurstStart = swrBurstInterval(s,1);
        Os(s).swrBurstEnd = swrBurstInterval(s,2);
        %     Os(s).starttimecheck0 = swrStart(s) - Os(s).starttime;
    end
end
%% add column: consuming
% fuck how do i do this? i need to keep the start time of the eclosing
% lickburst for each swr, such that i can use that start time to lookup
% proximity to a dio output time... i.e. i want to only keep the lickbouts
% that occurred close to when he received reward..
% alternative way is to timefilter with behavior indicator of correct state or
% not, but i'd have to check to make sure that the boundary between
% flipping between correct/incorrect is appropriate to capture the right
% bursts.. idk... i think maybe either would work, or do both and use to
% cross valididate
% ok now i can use the burst start/end + day, epoch of each ripple to
% lookup proximity to closest dioreward trigger.. 

dio = loaddatastruct(andef{2}, animal, 'DIO'); 
% outputdios = [7 8 9];
for s = 1:length(Os)
    dio{Os(s).day}{Os(s).epoch}
   
    
end

% need to definately check this with a data chunker.. 

% once i label each swr as being within a consumption lickbout or not, 
% determine if there's a difference in the swrlickphase distribution for
% the two

% then i can collect the other variables and do the testing of phase vs X.. corr

% and that will for now conclude the -behavioral testing of swr locking.. 


% then i can move on to doing the testing of burst vs nonburst
% swr=coactivity/replay score

% then i need to push all the animals through 




























