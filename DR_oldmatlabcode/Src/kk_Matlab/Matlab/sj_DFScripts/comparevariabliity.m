%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine','Ten'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
animals = {'Frank'};
%animals = {'Miles'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
%for i = 1:14
    
    %epochfilter{i} = ['(isequal($type,''sleep'')) & ($sleepnum == 2) & ($runbefore.exposureday == ',num2str(i),') & ($runbefore.dailyexposure == 1)'];  
%    epochfilter{i} = ['($exposureday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
    %epochfilter{i} = ['($exposureday == ',num2str(i),') & ($dailyexposure == 1)'];
    %epochfilter{i} = ['(($exposureday == ',num2str(i),') & ($dailyexposure == 1) & isequal($environment, ''TrackA''))|(($exposureday == ',num2str(i),') & ($exposureday < 9) & ($dailyexposure == 1) & isequal($environment, ''TrackB''))'];
%end

epochfilter{1} =  ['($exposureday == 1) & ($dailyexposure < 3) & isequal($environment, ''TrackA'')'];
epochfilter{2} =  ['($exposureday == 1) & ($dailyexposure < 3) & isequal($environment, ''TrackB'')'];

cellfilter = '((isequal($area, ''CA1'') | (isequal($area, ''CA3''))) && ($meanrate < 7) && ($numspikes >= 100))';

timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3) & ($mindisttowell >= 10))', 1},{'getriptimes','($nripples == 0)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))','minthresh',3} };



%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) < 3))', 6}, {'getriptimes','($nripples > 2)', [], 'cellfilter', '((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'}};
%timefilter = {{'getriptimes', '($nripples >= 1)', [], 'cellfilter', %'(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minthresh',3}};

%timefilter = {};

cf = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------




%run function- all cells within each day
%--------------------------------------------
iterator = 'multicellanal';

cf = setfilteriterator(cf,iterator);

cf = setfilterfunction(cf, 'calcplacevar', {'spikes', 'linpos', 'cellinfo'}, 5, 100);

cf = runfilter(cf);

pthresh = 0.01;
for a = 1:length(cf)
    for e = 1:length(epochfilter)
	for i = 1:length(cf(a).output{e})
	    o = cf(a).output{e}(i)
%	    sprintf('%s %d %d %d %d', cf(a).animal{1}, o.lf{1}.index);
	    ca1 = strcmp(o.cellloc, 'CA1');
	    tmp = max(find(ca1));
	    o.corr(find(o.corrp > pthresh)) = 0;
	    o.corr(find(~isfinite(o.corrp))) = 0;
	    imagesc(o.corr)
	    colorbar
	    hold on
	    if ~isempty(tmp)
		line([0 length(ca1)], [tmp+.5 tmp+.5]);
		line([tmp+.5 tmp+.5], [0 length(ca1)]);
	    end
	    pause
	    clf
	end
    end
end


