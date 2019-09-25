%Animal selection
%-----------------------------------------------------
animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};

%animals = {'Miles','Nine','Ten'};
%animals = {'Conley','Bond','Frank'};


%animals = {'Bond'};
%animals = {'Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------
epochfilters1 = [];
epochfilters2 = [];

epochfilter = [];
%epochfilter{1} = ['isequal($type,''sleep'')|isequal($type,''run'')'];

epochfilter1{1} = ['((isequal($type,''run'')) & ($dailyexposure == 1) & ($exposureday >6) & ($exposureday >6))'];  
epochfilter2{1} = ['((isequal($type,''sleep'')) & ($sleepnum == 2) & ($runbefore.exposureday >6) & ($runbefore.exposureday >6)  & ($runbefore.dailyexposure == 1))'];  


CA1f1 = [];
CA1f2 = [];
cmp = [];



%Filter creation
%--------------------------------------------------------

%CA1cellfilter = '($meanrate < 7)';
CA1cellfilter = '((isequal($area, ''CA1'')) && ($meanrate < 7))';
%CA3cellfilter = '((isequal($area, ''CA3'')) && ($meanrate < 7))';


%timefilter1 = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3} };
%timefilter2 = {{'getriptimes', '($nripples > 2)', [], 'cellfilter', '(isequal($area, ''CA1''))','minstd',3}};

timefilter1 = {{'get2dstate', '((abs($velocity) > 3))'},{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',2} };
%timefilter2 = {{'get2dstate', '((abs($velocity) < 1))'},{'getriptimes', '($nripples > 2)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};
%timefilter2 = {{'get2dstate', '((abs($velocity) < 1))'},{'getriptimes', '($nripples > 2)', [], 'cellfilter', '(isequal($area, ''CA1''))','minstd',3}};
%timefilter2 = {{'get2dstate', '(abs($velocity) < 1)'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};

timefilter2 = {{'get2dstate', '($immobilitytime > -1)'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};

%timefilter2 = {{'get2dstate', '($immobilitytime < 2.4) & (abs($velocity) < 1)'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};


CA1f1 = createfilter('animal',animals,'epochs',epochfilter1,'cells',CA1cellfilter,'excludetime', timefilter1);
CA1f2 = createfilter('animal',animals,'epochs',epochfilter2,'cells',CA1cellfilter,'excludetime', timefilter2);

%CA3f1 = createfilter('animal',animals,'epochs',epochfilter,'cells',CA3cellfilter,'excludetime', timefilter1);
%CA3f2 = createfilter('animal',animals,'epochs',epochfilter,'cells',CA3cellfilter,'excludetime', timefilter2);

%-----------------------------------------------------------

%Pairwise comparison
%-----------------------------------------------------------
iterator = 'singlecellanal';

CA1f1 = setfilteriterator(CA1f1,iterator);
CA1f2 = setfilteriterator(CA1f2,iterator);
%CA3f1 = setfilteriterator(CA3f1,iterator);
%CA3f2 = setfilteriterator(CA3f2,iterator);


%CA1f1 = setfilterfunction(CA1f1, 'calcpeakrate', {'spikes', 'linpos'});
CA1f1 = setfilterfunction(CA1f1, 'calctotalmeanrate', {'spikes'},'appendindex',1);
CA1f2 = setfilterfunction(CA1f2, 'calctotalmeanrate', {'spikes'},'appendindex',1);

%CA3f1 = setfilterfunction(CA3f1, 'calctotalmeanrate', {'spikes'},'appendindex',1);
%CA3f2 = setfilterfunction(CA3f2, 'calctotalmeanrate', {'spikes'},'appendindex',1);


%CA1f1 = setfilterfunction(CA1f1, 'calcrippleactivationprob', {'spikes','ripples','cellinfo'},'appendindex',1,'ratio',1);
%CA1f2 = setfilterfunction(CA1f2, 'calcrippleactivationprob', {'spikes','ripples','cellinfo'},'appendindex',1,'ratio',1);

CA1f1 = runfilter(CA1f1);
CA1f2 = runfilter(CA1f2);
%CA3f1 = runfilter(CA3f1);
%CA3f2 = runfilter(CA3f2);

CA1f1groups = numericgroupcombine(CA1f1,1);
CA1f2groups = numericgroupcombine(CA1f2,1);
%CA3f1groups = numericgroupcombine(CA3f1,1);
%CA3f2groups = numericgroupcombine(CA3f2,1);

%indexcolumns = [1 2 3 4 5];

CA1f1groups{1} = CA1f1groups{1}(find(CA1f1groups{1}(:,3)==2),:)
indexcolumns = [1 2 4 5];

cmp{1} = indexmatch(CA1f1groups{1},CA1f2groups{1},indexcolumns);
%cmp{1} = indexmatch(CA3f1groups{1},CA3f2groups{1},indexcolumns);

plot(cmp{1}(:,5),cmp{1}(:,6),'.')
[a,b] = corrcoef(cmp{1}(:,5),cmp{1}(:,6));

%--------------------------------------------------------------
%2-contigencies

goodcells1 = cmp{1}((cmp{1}(:,3) == 2),:);
goodcells2 = cmp{1}(((cmp{1}(:,3) == 3)&(cmp{1}(:,6) < .1)),:);
%goodcells2 = cmp{1}(((cmp{1}(:,3) == 3)),:);
indexcolumns2 = [1 2 4 5];
result = indexmatch(goodcells1,goodcells2,indexcolumns2);


%column 5,7 is no ripples
%column 6,8 is during ripples
plotcolumns = [5 8]

X = result(:,plotcolumns(1));
%X = xdata;
X(:,2) = 1;
Y = result(:,plotcolumns(2));
[b,bint,r,rint,stats] = regress(Y,X);
stats
 figure
 plot(result(:,plotcolumns(1)),result(:,plotcolumns(2)),'.')
 
%-----------------------------------------------------------

%3-contingencies
% goodcells1 = cmp{1}((cmp{1}(:,3) == 6),:);
% goodcells2 = cmp{1}(((cmp{1}(:,3) == 2)&(cmp{1}(:,6) < .1)),1:5);
% indexcolumns = [1 2 4 5];
% result = indexmatch(goodcells1,goodcells2,indexcolumns);
% goodcells3 = cmp{1}(((cmp{1}(:,3) == 7)&(cmp{1}(:,6) < .1)),[1 2 4:end]);
% indexcolumns = [1 2 3 4];
% result = indexmatch(result,goodcells3,indexcolumns);
% 

% X = result(:,6);
% X(:,2) = 1;
% Y = result(:,8);
% [b,bint,r,rint,stats] = regress(Y,X);

%  figure
%  plot(result(:,6),result(:,8),'.')

 

