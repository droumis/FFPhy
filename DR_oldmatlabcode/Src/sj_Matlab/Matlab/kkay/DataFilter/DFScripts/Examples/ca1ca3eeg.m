
%Animal selection
%-----------------------------------------------------
animals = {'Frank'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
% only examine the first exposure to a novel environment
epochfilter{1} = ['($exposure == 1)'];

% select pairs of tetrodes in CA3
tetrodepairfilter = {'(isequal($area, ''CA3''))', '(isequal($area, ''CA3''))'};

% exclude ripples and low velocity times
timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}};


iterator = 'eeganal';

f = createfilter('animal',animals,'epochs',epochfilter,'tetrodepairs',tetrodepairfilter,'excludetimefilter', timefilter, 'iterator', iterator);

% calculate the phase difference between the theta waveforms
f = setfilterfunction(f, 'calcphasediff', {'theta'});
f = runfilter(f);
