
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
epochfilter{1} = ['($exposure == 6)'];
%epochfilter{2} = ['($exposure == 6)'];

%cellpairfilter = {'allcomb', '(isequal($area, ''CA1'') && ($meanrate < 7))', '(isequal($area, ''CA3'') && ($meanrate < 7))'};
cellpairfilter = {'allcomb', '(isequal($area, ''CA3'') && ($meanrate < 7))', '(isequal($area, ''CA1'') && ($meanrate < 7))'};

%timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6}};
timefilter = {{'gethighthetatimes', '($nhightheta > 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};


iterator = 'singlecellanal';

f = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcxcorrmeasures', {'spikes', 'linpos'}, 'edgespikes', 1, 'calclinfields', 1, 'calctrajxcorr', 1);
f = runfilter(f);

% plot correlations and fields
g = gaussian(3, 18);
phbins = [0:pi/12:2*pi];
for a = 1:length(f)
    for e = 1:length(f(a).output)
	for c = 1:length(f(a).output{e})
	    xc = f(a).output{e}(c);
	    if (~isempty(xc.c1vsc2) & (max(xc.ac1) > 100) & (max(xc.ac2) > 100))
		figure(1);
		for i = 1:4
		    subplot(4,1,i);
		    tmpind = find(abs(xc.time) < .5);
		    if (~isempty(xc.c1vsc2{i}))
			plot(xc.time(tmpind), smoothvect(xc.c1vsc2{i}(tmpind), g));
		    end
		    if (i == 1)
			title(sprintf('%d %d - %d %d vs. %d %d', xc.index))
		    end
		end
		figure(2);
		lf1 = xc.lf1.trajdata;
		lf2 = xc.lf2.trajdata;
		pp1 = xc.lf1.phasedist;
		pp2 = xc.lf2.phasedist;
		ph1 = [];
		ph2 = [];
		for i = 1:4
		    subplot(4,2,2*i-1);
		    plot(lf1{i}(:,1), lf1{i}(:,5), 'b');
		    hold on;
		    plot(lf2{i}(:,1), lf2{i}(:,5), 'r');
		    set(gca, 'XLim', [0 180]);
		    subplot(4,2,2*i);
		    plot(pp1{i}.dist, pp1{i}.phase, 'b.');
		    hold on;
		    plot(pp1{i}.dist, pp1{i}.phase+2*pi, 'b.');
		    plot(pp2{i}.dist, pp2{i}.phase, 'r.');
		    plot(pp2{i}.dist, pp2{i}.phase+2*pi, 'r.');
		    ph1 = [ph1 ;pp1{i}.phase];
		    ph2 = [ph2 ;pp2{i}.phase];
		    set(gca, 'YLim', [0 4*pi]);
		    set(gca, 'XLim', [0 180]);
		end
		figure(3);
		subplot(2,1,1)
		c = histc(ph1, phbins);
		bar(phbins, c, 'b');
		subplot(2,1,2)
		c = histc(ph2, phbins);
		bar(phbins, c, 'r');
		'pausing'
		pause
		clf(1);
		clf(2);
	    end
	end
    end
end
	




