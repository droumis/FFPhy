
%{

- need to finish XP mod SU figure
.. what is the actual final result?
.. where is the result saved? is there an accessible record of per cell
results?
.. what is the final result of the swr mod su?

-- behavior position fig.. 

%}
loaddata = 0;
pconf = paramconfig;

create_filter = 0;
run_ff = 0;
load_ffdata = 1;

make_expvarCat = 0;
load_expvarCat = 0;

calcSUphasemod = 
loadSUPhaseMod = 
gatherPhaseModResults = 

% data filter params
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_eventTrigSpiking';
expvars = {'all', 'wetLickBursts', 'dryLickBursts'};

Fp.Label = 'wtrackLickTrigSpiking';
Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
    'excludeAfterLastWell', 'nonMU_cells', Fp.Label, Fp.filtfunction};

%% FF
Fp = load_filter_params(Fp);
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes',...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator, 'cells',...
        Fp.cellfilter);
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

%% make design mat to slice trials
if make_expvarCat
    dmat = makeExpvarCatDesignMat(F, expvars, 'eventType', Fp.eventType);
end

if load_expvarCat
    outdir = 'expvarCat';
    outpath = [pconf.andef{2},outdir,'/'];
    dmat = load_data(outpath, [outdir,'_',Fp.env,'_',Fp.eventType], Fp.animals);
end

%% Spike Phase mod
if calcSUphasemod
    pmodF = calcPhaseMod(F, dmat); % Saw
    save_data(pmodF, 'results', [Fp.Label '_phasemod']);
end
if loadSUPhaseMod
    pmodF = load_data('results', [Fp.Label '_phasemod'], Fp.animals);
end

%% SPIKE Gather all animals Phasemod su per area, eventSet
if gatherPhaseModResults
    areaPhasemod = {};
    allAnPmodF = cell2mat(arrayfun(@(x) pmodF(x).output{1}', 1:length(pmodF), 'un', 0)');
    numESet = length(pmodF(1).dmatIdx);
    
    for ar = 1:length(Fp.areas) % per area
        % find cells in this area
        areaIdx = strcmp({allAnPmodF.area}', Fp.areas{ar}{1});
        subareaIdx = ~cellfun(@isempty, strfind({allAnPmodF.subarea}', Fp.areas{ar}{2}), 'un', 1);
        iareaIdx = find(all([areaIdx subareaIdx],2));
        iareaIdx = iareaIdx(arrayfun(@(x) ~isempty(allAnPmodF(x).phasemod), ...
            iareaIdx,'un',1));
        
        for iv = 1:numESet % per eventSet
            gud = [];
            dumpy = [];
            phasemod = [];
            mPctChangeSh = [];
            for i = 1:length(iareaIdx)
                try
                    m = allAnPmodF(iareaIdx(i)).phasemod{iv};
                    %                     mSh = allAnPmodF(iareaIdx(i)).mPctChangeSh{iv};
                    if ~isempty(m)
                        phasemod = [phasemod; m];
                        mPctChangeSh = [mPctChangeSh; mSh];
                        gud = [gud; i];
                    end
                catch
                    continue
                end
            end
            areaPhasemod{ar,iv} = phasemod;
        end
    end
end


