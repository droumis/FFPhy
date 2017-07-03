
% get spikes and pos

% warning('off','all');
close all
runFilterFramework = 1;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;

savefigs= 0;
pausefigs = 0;
investInfo = animaldef(lower('Demetris'));
%% ---------------- plotting params --------------------------
colorSet = 'DR1';
fonttype = 'Arial';
titlesize = 16;
axissize = 16;
arrowsize = 12;
usecolormap = 'hot';
trajname = {'outbound A', 'inbound A', 'outbound B', 'inbound B'};
SpacingHoriz = 0.02;
SpacingVert = 0.02;
Padding = 0.00;
position = [.1 .1 .8 .5];
Spacing = 0.00;
Padding = 0.00;
MarginLeft = 0.05;
MarginRight = 0.05;
MarginTop = 0.15;
MarginBottom =  0.1;
%% ---------------- Data Filters --------------------------
animals = {'JZ1'};
days = [4];
epochType = {'run'};
epochEnvironment = {'openfield'}; %wtrack

ntAreas = {'ca1', 'mec', 'por', 'v2l', 'sub'};
filtfunction = 'spikePos';
eventtype = 'rippleskons';
eventSourceArea = 'ca1';
consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 5;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minripvelocity = 0;
maxripvelocity = 4;
run2D = 1;
% runLin = 1;
%% ---------------- Paths and Title strings ---------------------------------------------------
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
resultsOutDirectory = sprintf('%s%s/', investInfo{3}, filtfunction);
figdirectory = sprintf('%s%s/', investInfo{4}, filtfunction);
filenamesave = sprintf('%s_D%s', strjoin(epochEnvironment,'-'), strjoin(arrayfun(@(x) num2str(x),days,'UniformOutput',false),'-')); %add more naming stuff
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
filenameTitle = strrep(filename,'_', '\_');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1
    eptypeEnv = [epochType; epochEnvironment];
    epochfilter = sprintf('((isequal($type, ''%s'')) && (isequal($environment, ''%s''))) || ',eptypeEnv{:});
    yuck = strfind(epochfilter, '||');
    epochfilter = epochfilter(1:yuck(end)-1);  %cut off the trailing '||'
    cellfilter = '($numspikes > 10)';
    tetfilter = sprintf('(isequal($area,''%s'')) || ', ntAreas{:});
    yuck = strfind(tetfilter, '||');
    tetfilter = tetfilter(1:yuck(end)-1); %cut off the trailing '||'
    timefilter{1} = {'get2dstate', '(abs($velocity) >= 3)'};
    timefilter{2} = {'getconstimes', '($cons == 0)', [eventSourceArea eventtype],1,'consensus_numtets',consensus_numtets,...
        'minstdthresh',minstdthresh,'exclusion_dur',exclusion_dur,'minvelocity',minripvelocity,'maxvelocity',maxripvelocity};
    iterator = 'singlecellanal';
    if run2D
        F = createfilter('animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
        F = setfilterfunction(F, ['dfa_' filtfunction], {'spikes', 'linpos', 'pos', 'task'});
        F = runfilter(F);
    end
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(filtOutputDirectory);
        mkdir(filtOutputDirectory);
    end
    save(sprintf('%s/%s',filtOutputDirectory, filename), 'F', '-v7.3');
    disp(sprintf('filteroutput saved to %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/%s',filtOutputDirectory, filename))
    disp(sprintf('filteroutput loaded: %s/%s',filtOutputDirectory, filename))
end


if ~isdir(resultsOutDirectory);
        mkdir(resultsOutDirectory);
    end
save(sprintf('%s/%s',resultsOutDirectory, filename), 'JZ1spikesPosD4', '-v7.3');
disp(sprintf('results out saved to: %s/%s',resultsOutDirectory, filename));
load(sprintf('%s/%s',resultsOutDirectory, filename));