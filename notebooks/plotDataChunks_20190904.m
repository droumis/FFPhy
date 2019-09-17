

%{
Change this to either plot all datachunks centered on each ripple..
OR plot all the 10s timechunks, regaredless of whether there were ripples
during it or not
%}
clear all
create_filter = 1;
run_ff = 1;
save_ffdata = 0;
load_ffdata = 0;

loadEpoch = 0;
pilotEpoch = 0;

%% data filter params
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_plotDataChunks';
Fp.add_params = {'savefigs', 'wtrackdays','excludePriorFirstWell','<4cm/s', ...
    'excludeAfterLastWell', 'valid_ntrodes', 'exemplar_wepochs'}; %'excludeNoise',
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
%% FF
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator,'eegtetrodes',Fp.tetfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F); for d = 1:length(F); F(d).datafilter_params = Fp; end; end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s', Fp.epochEnvironment)); end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s', Fp.epochEnvironment)); end

%% pilot one epoch.. bypasses iterator
if pilotEpoch
    if loadEpoch
        day = F(1).epochs{1}(1,1);
        epoch = F(1).epochs{1}(1,2);
        % replaces the iterator to load a single epoch
        animal = Fp.animals{1};
        animdef = animaldef(animal);
        ntrodes = F.eegdata{1}{1};
        events = loaddatastruct(animdef{2}, animal, 'ca1rippleskons', day);
        %     noise = loaddatastruct(animdef{2}, animal, 'ca1noisekons', day);
        eeg = loadeegstruct(animdef{2}, animal, 'eeg', day, epoch, ntrodes);
        ripple = loadeegstruct(animdef{2}, animal, 'ripple', day, epoch, ntrodes);
        theta = loadeegstruct(animdef{2}, animal, 'theta', day, epoch, ntrodes);
        pos = loaddatastruct(animdef{2}, animal, 'pos', day);
        linpos = loaddatastruct(animdef{2}, animal, 'linpos', day);
        spikes = loaddatastruct(animdef{2}, animal, 'spikes', day);
        cellinfo = loaddatastruct(animdef{2}, animal, 'cellinfo', day);
        tetinfo = loaddatastruct(animdef{2}, animal, 'tetinfo', day);
        DIO = loaddatastruct(animdef{2}, animal, 'DIO', day);
        task = loaddatastruct(animdef{2}, animal, 'task', day);
    end
    minstdthresh = Fp.minstdthresh;
    maxvelocity = Fp.maxvelocity;
    epidx = find(ismember(F.epochs{1}, [day epoch], 'rows'));
    excludeIntervals = F.excludetime{1}{epidx};
    dfa_plotDataChunks([day epoch], excludeIntervals, 'ca1rippleskons', events, ...
        'eeg', eeg, 'ripple', ripple, 'theta', theta, 'pos', pos, 'linpos', linpos, ...
        'spikes', spikes, 'cellinfo', cellinfo, 'tetinfo', tetinfo, 'DIO', DIO, ...
        'task', task, 'minstdthresh', Fp.minstdthresh, 'maxvelocity', Fp.maxvelocity, ...
        'pausefigs', pausefigs, 'savefigs', savefigs, 'centerEvents', Fp.centerEvents, ...
        'splitSize', Fp.splitSize, 'animal', animal, 'skipExist', Fp.skipExist, 'Yoffset', ...
        Fp.Yoffset, 'savefigas', Fp.savefigas);
end
