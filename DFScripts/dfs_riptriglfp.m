
% clear all
close all
runFilterFramework = 0;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 1;
% EpochMean = 1;
resaveFilterOutput = 1;
% plotSpec_EpochMean = 1;
% plotSpec_allEpochs = 0;
outputDirectory = '/typhoon/droumis/analysis';
%% ---------------- plotting params --------------------------


%% ---------------- Data Filters --------------------------

animals = {'D13'};
days = [1];
eventtype = 'rippleskons';
% eventarea = 'ca1';
epochEnvironment = 'wtrack';% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
epochType = 'run';
% tetAreas = ['ca1', 'mec', 'por']; %ca1, mec, por

consensus_numtets = 1;   % minimum # of tets for consensus event detection
minthresh = 2;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;

filename = sprintf('riptriglfp_%s%s_%s_%sTets_%s.mat', eventarea, eventtype, epochEnvironment, tetAreas, cell2mat(animals));
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
    epochfilter =    sprintf('(isequal($type, ''%s'')) && (isequal($environment, ''%s''))',epochType, epochEnvironment); %'isequal($type, ''run'') && (isequal($environment, ''MultipleW''))'; %%'(isequal($type, ''sleep''))'; %%%&& isequal($descript, ''post-allruns''))';%   %%% %'isequal($type, ''run'') && isequal($environment, ''WTrackA'') && (($exposure>=1) && ($exposure<=10))';  %
    iterator = 'epocheeganal';
    tetfilter = '(isequal($area,''ca1'') || isequal($area,''mec'') || isequal($area,''por''))'; % || isequal($area,''v2l'') || isequal($area,''sub''))';
    % timefilter{1} = {'get2dstate','($velocity<4)'};
    % timefilter{2} = {'kk_getriptimes','($nripples>=1)',[],'tetfilter',tetfilter,'minthresh',5};
    timefilter{1} = {'getconstimes', '($cons == 1)',[eventarea,eventtype],1,'consensus_numtets',consensus_numtets,...
        'minthresh',minthresh,'exclusion_dur',exclusion_dur,'minvelocity',minvelocity,'maxvelocity',maxvelocity};
    %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
    F = createfilter('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
    %----------f = setfilteriterator(f, funcname, loadvariabl   es, options)--------
    F = setfilterfunction(F, 'dfa_riptriglfp', {'eeg', [eventarea,eventtype]},'eventtype',eventtype);
    tic
    F = runfilter(F);
    F.filterTimer = toc;
    F.worldDateTime = clock;
    F.dataFilter = struct('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
    
    %% ---------------- Collect across tetrodes ---------------------------------------------------
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(sprintf('%s/filter_output/riptrigspectrogram/', outputDirectory));
        mkdir(sprintf('%s/filter_output/riptrigspectrogram/', outputDirectory));
    end
    %save the entire workspace for filter provenance
    save(sprintf('%s/filter_output/riptrigspectrogram/riptrigspec_%s%s_%s_%sTets_%s.mat',outputDirectory, eventarea, eventtype, epochEnvironment, tetAreas, cell2mat(animals)), 'F','-v7.3');
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
%     load(sprintf('%s/filter_output/riptrigspectrogram/riptrigspec_%s%s_%s_%sTets_%s.mat',outputDirectory, eventarea, eventtype, epochEnvironment, tetArea, cell2mat(animals)));
    load(sprintf('%s/filter_output/riptrigspectrogram/%s',outputDirectory,filename))
end