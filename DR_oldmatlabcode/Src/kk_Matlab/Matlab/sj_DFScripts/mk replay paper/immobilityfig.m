%Animal selection
%-----------------------------------------------------
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};

%animals = {'Miles','Nine','Ten'};
animals = {'Conley','Bond','Frank'};
%animals = {'Frank'};

%animals = {'Bond'};
%animals = {'Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------

epochfilter = [];
epochfilter{1} = ['isequal($type,''sleep'')'];
%epochfilter{1} = ['isequal($type,''sleep'')|isequal($type,''run'')'];
%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];  


f1 = [];
f2 = [];
cmp = [];



%Filter creation
%--------------------------------------------------------

cellpairfilter = {'allcomb', '(isequal($area, ''CA1'') && ($meanrate < 7))', '(isequal($area, ''CA1'') && ($meanrate < 7))'};
%cellpairfilter = {'allcomb', '(($meanrate < 7))', '(($meanrate < 7))'};


%timefilter1 = {{'get2dstate', '((abs($velocity) > 1))'},{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',2} };
%timefilter2 = {{'get2dstate', '((abs($velocity) < 1))'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};

timefilter1 = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minthresh',2} };
%timefilter2 = {{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};

timefilter2 = {{'get2dstate', '($immobilitytime == 0)'}};
%timefilter2 = {{'gethighthetatimes', '($nhightheta == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};
%timefilter2 = {};


%timefilter2 = {{'get2dstate', '((abs($velocity) < 3))'}};
%timefilter2 = {{'getlinstate', '($distance2well < 1.5)',6},{'getriptimes', '($nripples > 2)', [], 'cellfilter', '(isequal($area, ''CA1''))','minstd',3}};

f2 = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetime', timefilter2);

% go through all of the time filters and collect the length of the exclude
% windows
exlen = [];
for a = 1:length(f2)
    for e = 1:length(f2(a).excludetime)
	for p = 1:length(f2(a).excludetime{e})
	    exlen = [exlen ; f2(a).excludetime{e}{p}(:,2) -  ...
	                     f2(a).excludetime{e}{p}(:,1)];
	end
    end
end

cutoff = 5;
subplot(2,1,1)
[a b] = flhist(exlen(find(exlen < cutoff)), 20);
bar(b,a, 'FaceColor', [0 0 0]);

subplot(2,1,2)
bins = 5:20:max(exlen);
[a b] = flhist(exlen(find(exlen >= cutoff)), bins);
bar(b,a, 'FaceColor', [0 0 0]);
