
run('setparameters');

%-----------------------------------------------------------
fprintf('Generating filters...');
%-----------------------------------------------------------
pulsefilter = {'(round($amplitude) == 200)'};

tetrodes = [1:32];
tetrodefilter = '$chan>=0';

% (FILTER ORDER IS IMPORTANT B/C OF USE OF SHORT CIRCUIT && !!
% epochfilter = 'strcmp($type,''run'') && ($exposure(1)==1) && ($exposure(2)>=1)';
% epochfilter = '($epoch==9)';
% epochfilter = 'strcmp($type,''run'') && ($day==3) && ($epoch==8)';
epochfilter = 'strcmp($type,''run'')';

runEpochs = createfilter('animal',animals,'epochs',epochfilter);
runEpochs = combinefilterepochs(runEpochs,'split','empty'); % split all multiepochs clear empty epochs
% runEpochs = combinefilterepochs(runEpochs,'days','empty'); % combine same days and clear empty epochs


runEEG = modifyfilter(runEpochs, 'tetrodes', tetrodefilter, ...
  'iterator', 'timewindowanal', 'stim', pulsefilter, ...
  'filterfunction', {'getraweegwindow', {'eegfiledir'}, 'eegwindow', pulsewin});

runBehav = modifyfilter(runEEG, 'iterator', 'epochbehavtimeanal', ...
  'filterfunction', {'getbehavstatetime',{'pos','linpos','task','track_regions'}, ...
                     'behavVars',{'velocity','position'}});

%-----------------------------------------------------------
fprintf('\nRunning filters...');
%-----------------------------------------------------
runEEG = runfilter(runEEG);
runBehav = runfilter(runBehav);
fprintf('\n');

