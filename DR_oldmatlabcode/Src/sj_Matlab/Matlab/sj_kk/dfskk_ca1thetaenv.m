
%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Ten'};
%animals = {'Dudley'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
animals = {'Frank'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = '(isequal($type, ''run''))';
%epochfilter{2} = ['($exposure == 6)'];

ca1tetfilter = '(isequal($area, ''CA1''))';
ca3tetfilter = '(isequal($area, ''CA3''))';


timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}};


iterator = 'eeganal';

f = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',ca3tetfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcenvdist', {'theta'});
f = runfilter(f);

'pause'

% plot the histograms
for a = 1:length(f)
    for e = 1:length(f(a).output)
	for t = 1:length(f(a).output{e})
	    bar(f(a).output{e}(t).bins, f(a).output{e}(t).envhist) 
	    pause
	    clf
	end
    end
end
