
%Animal selection
%-----------------------------------------------------
animals = {'Bond','Frank','Ten'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter{1} = ['isequal($environment, ''TrackA'') | isequal($environment, ''TrackB'')'];
%for i = 1:14
    %epochfilter{i} = ['isequal($environment, ''TrackA'') | isequal($environment, ''TrackB'')'];
%   epochfilter{i} = ['isequal($environment, ''TrackA'')'];
%nd

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7) && ($propbursts < 0.1) && ($csi < .05) && ($numspikes > 200))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7) && ($propbursts < 0.1) && ($csi < .05) && ($numspikes > 200))';

%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 0))', 6} };
timefilter = {{'getriptimes', '($nripples > 1)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};

%timefilter = {};

iterator = 'singlecellanal';
%iterator = 'multicellanal';

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter, 'iterator', iterator);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter, 'iterator', iterator);
%-----------------------------------------------------------




%ca1f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes', 'linpos'});
ca1f = setfilterfunction(ca1f, 'calctotalmeanrate', {'spikes'});
%ca3f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes', 'linpos'});

ca1f = runfilter(ca1f);
%ca3f = runfilter(ca3f);
 
tmp = numericgroupcombine(ca1f,1);
'pause'
pause

% plot the ca1 ivy cells
for a = 1:length(ca1f)
    if (~isempty(ca1f(a).output))
	for i = 1:length(ca1f(a).output{1})
	    trajdata = ca1f(a).output{1}(i).trajdata;
	    for t = 1:4
		subplot(4,1,t)
		plot(trajdata{t}(:,1), trajdata{t}(:,5));
	    end
	    title(sprintf('%d %d %d %d', ca1f(a).output{1}(i).index));
	    pause
	end
    end
end


for a = 1:length(ca3f)
    if (~isempty(ca3f(a).output))
	for i = 1:length(ca3f(a).output{1})
	    trajdata = ca3f(a).output{1}(i).trajdata;
	    for t = 1:4
		subplot(4,1,t)
		plot(trajdata{t}(:,1), trajdata{t}(:,5));
	    end
	    title(sprintf('%d %d %d %d', ca3f(a).output{1}(i).index));
	    pause
	end
    end
end


