
%{
 plot lfp, rips, pos, speed

- label/color ntrodes, organize by area/layer
- filter ripples for only the included ones

maybe i should split all the data into 100 second intervals? 
i think i've done this before, no?
- save out as matlab figures so i can load and zoom in on each minute
%}

savefigs = 0;
pausefigs = 1;

% animal = 'D10';
% animdef = animaldef(animal);
% day = 4;
% epoch = 2;
% splitSize = 10; % seconds
% Yoffset = 600; % lfp y plot offset

Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ4'}; %;{'D10'}; %
% Fp.animals = {'JZ3'};
Fp.filtfunction = 'dfa_plotDataChunks';
Fp.add_params = {'wtrackdays','excludeNoise','excludePriorFirstWell','<4cm/s', ...
    'excludeAfterLastWell', 'nonref_ntrodes'};
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
pconf = paramconfig;
%%
make_swrLFP = 1;
save_swrLFP = make_swrLFP;
load_swrLFP = 0;
if make_swrLFP
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter, 'iterator', Fp.iterator, 'eegtetrodes', Fp.tetfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    F = runfilter(F); for d = 1:length(F); F(d).datafilter_params = Fp; end; end
if save_swrLFP
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s', Fp.epochEnvironment)); end
if load_swrLFP
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s', Fp.epochEnvironment)); end

% for ade = 1:length([1])

%% load data
lfp = loadeegstruct(animdef{2}, animal, 'eeg', day, epoch, 1:30);
swr = loaddatastruct(animdef{2}, animal, 'ca1rippleskons');
pos = loaddatastruct(animdef{2}, animal, 'pos');
