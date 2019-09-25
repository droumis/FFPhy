
run('setparameters');

%-----------------------------------------------------------
fprintf('Generating filters...');
%-----------------------------------------------------------
pulsefilter = {'(round($timeuntilnext) >= 10000) & (round($timesincelast) >= 10000) & (round($pulselength/10) == min(unique(round($pulselength/10))))'};

tetrodes = [1:32];
tetrodefilter = '$tetrode>=0';

% (FILTER ORDER IS IMPORTANT B/C OF USE OF SHORT CIRCUIT && !!
epochfilter = 'strcmp($type,''run'') && ($exposure(1)==1) && ($exposure(2)>=1)';

runEpochs = createfilter('animal',animals,'epochs',epochfilter);
runEpochs = combinefilterepochs(runEpochs,'split','empty'); % split all multiepochs clear empty epochs
% runEpochs = combinefilterepochs(runEpochs,'days','empty'); % combine same days and clear empty epochs


runEEG = modifyfilter(runEpochs, 'tetrodes', tetrodefilter, ...
  'iterator', 'timewindowanal', 'stim', pulsefilter, ...
  'filterfunction', {'getraweegwindow', {'eegfiledir'}, 'eegwindow', pulsewin});

runBehav = modifyfilter(runEEG, 'iterator', 'epochbehavtimeanal', ...
  'filterfunction', {'getbehavstatetime',{'pos','linpos','task','track_regions'}, ...
                     'behavVars',{'smoothvelocity','position'},'t_offset',0});

runVelWindow = modifyfilter(runEEG, 'iterator', 'epochbehavtimeanal', ...
  'filterfunction', {'getbehavstatewindow',{'pos','linpos','task','track_regions'}, ...
                     'behavVar','smoothvelocity','behavwindow',speedwin});

runMulti = modifyfilter(runEpochs, 'tetrodes', tetrodefilter, ...
  'iterator', 'timewindowanal', 'stim', pulsefilter, ...
  'filterfunction', {'getspikewindow', {'multi'}, ...
    'spikewin', spikewin, 'output','psth', ...
    'inputType','multi','smoothing', ''});

runTheta = modifyfilter(runEpochs, 'tetrodes', tetrodefilter, ...
  'iterator', 'timewindowanal', 'stim', pulsefilter, ...
  'filterfunction', {'geteegwindow', {'theta'}, 'eegwindow', pulsewin});


%-----------------------------------------------------------
fprintf('\nRunning filters...');
%-----------------------------------------------------
fprintf(' [E]'); runEEG = runfilter(runEEG);
fprintf(' [B]'); runBehav = runfilter(runBehav);
fprintf(' [V]'); runVelWindow = runfilter(runVelWindow);

fprintf(' [M]'); runMulti = runfilter(runMulti);
fprintf(' [T]'); runTheta = runfilter(runTheta);
fprintf('\n');

