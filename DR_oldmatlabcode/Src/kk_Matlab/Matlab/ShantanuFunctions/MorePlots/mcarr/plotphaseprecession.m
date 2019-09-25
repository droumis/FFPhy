%RUN PLOT PHASE PRESESSION
animals = {'Eight'};

% Epoch selection
epochfilter = [];
epochfilter{1} = ['$exposure == 10 & isequal($description,''TrackA'')'];


%Time selection
timefilter = {{'getriptimes', '($nripples == 0)', [], 'tetfilter', '(isequal($area, ''CA1''))'},{'getlinstate', '($traj == 1|3 & (abs($velocity) >= 3))', 6}};

% Cell selection
ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7) && $numspikes>200)';

%Select iterator
iterator = 'singlecelleeganal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'cells',ca1cellfilter,'iterator', iterator);
f = setfilterfunction(f, 'plotphaseprecession', {'spikes','linpos','eeg'});
f = runfilter(f);
