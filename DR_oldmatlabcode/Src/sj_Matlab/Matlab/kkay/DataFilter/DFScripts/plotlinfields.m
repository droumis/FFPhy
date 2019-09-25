
%Animal selection
%-----------------------------------------------------
%animals = {'Bond','Frank','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Frank'};
animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
%epochfilter{1} = ['isequal($environment, ''TrackB'') && ($dailyexposure == 1)'];
epochfilter{1} = ['isequal($environment, ''TrackA'') && ($dailyexposure == 1)'];
%for i = 1:14
    %epochfilter{i} = ['isequal($environment, ''TrackA'') | isequal($environment, ''TrackB'')'];
%   epochfilter{i} = ['isequal($environment, ''TrackA'')'];
%nd

cellfilter = '($meanrate < 7)';
ca1cellfilter = '(isequal($area, ''CA1''))';
ca3cellfilter = '(isequal($area, ''CA3''))';

%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 0))', 6} };
timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};

%timefilter = {};

iterator = 'singlecellanal';
%iterator = 'multicellanal';

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'iterator', iterator);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'iterator', iterator);
%-----------------------------------------------------------
cf = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'iterator', iterator);
%-----------------------------------------------------------




%ca1f = setfilterfunction(ca1f, 'cellprop', {'spikes'});
%ca3f = setfilterfunction(ca3f, 'cellprop', {'spikes'});
ca3f = setfilterfunction(ca3f, 'filtercalclinfields', {'spikes', 'linpos'});
ca1f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes', 'linpos'});
cf = setfilterfunction(cf, 'filtercalclinfields', {'spikes', 'linpos'});

'ca1'
%ca1f = runfilter(ca1f);

'ca3'
%ca3f = runfilter(ca3f);
 
cf = runfilter(cf);

figure
for i = 1:length(cf.output{1})
    r = cf.output{1}(i);
    for s = 1:4
	subplot(4,1,s);
	plot(r.trajdata{s}(:,1), r.trajdata{s}(:,5));
    end
    subplot(4,1,1);
    title(sprintf('%d %d %d %d', r.index));
    pause
    clf
end


