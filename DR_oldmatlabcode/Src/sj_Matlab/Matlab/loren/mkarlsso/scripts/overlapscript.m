%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Ten','Alex','Bond','Frank','Nine'};
animals = {'Miles','Conley','Ten','Bond','Frank','Nine'};


%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackA'')'];
epochfilter{2} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];

%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];
%epochfilter{2} = ['($dailyexposure == 2) & isequal($environment, ''TrackB'')']; %use only the first day days in the novel track

ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 7))';

timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6} };
%timefilter = {};

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter);
%-----------------------------------------------------------

%Pairwise comparison
%-----------------------------------------------------------
iterator = 'singlecellanal';

ca1f = setfilteriterator(ca1f,iterator);
ca3f = setfilteriterator(ca3f,iterator);

ca1f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes','linpos'});
ca3f = setfilterfunction(ca3f, 'filtercalclinfields', {'spikes','linpos'});

ca1f = runfilter(ca1f);
ca3f = runfilter(ca3f);

ca1overlap = [];
for i = 1:length(ca1f)
    ca1overlap{i} = [];
    if (length(ca1f(i).output) > 1)
        g1index = reshape([ca1f(i).output{1}.index]',4,[])';
        g2index = reshape([ca1f(i).output{2}.index]',4,[])';
        matches = rowfind(g1index(:,[1 3 4]),g2index(:,[1 3 4]));
        %matches = matches(randperm(length(matches))); %random permutation
        for j = 1:length(matches)
            if (matches(j) > 0)
                tmpoverlap = calcoverlap(ca1f(i).output{1}(j).trajdata, ca1f(i).output{2}(matches(j)).trajdata,'thresh',3);
                if ~isempty(tmpoverlap)
                    ca1overlap{i} = [ca1overlap{i}; [i ca1f(i).output{1}(j).index([1 3 4]) tmpoverlap]];
                end
            end
        end
    end
end

ca3overlap = [];
for i = 1:length(ca3f)
    ca3overlap{i} = [];
    if (length(ca3f(i).output) > 1)
        g1index = reshape([ca3f(i).output{1}.index]',4,[])';
        g2index = reshape([ca3f(i).output{2}.index]',4,[])';
        matches = rowfind(g1index(:,[1 3 4]),g2index(:,[1 3 4]));
        %matches = matches(randperm(length(matches))); %random permutation
        for j = 1:length(matches)
            if (matches(j) > 0)
                tmpoverlap = calcoverlap(ca3f(i).output{1}(j).trajdata, ca3f(i).output{2}(matches(j)).trajdata,'thresh',3);          
                
                if ~isempty(tmpoverlap)
                    ca3overlap{i} = [ca3overlap{i}; [i ca3f(i).output{1}(j).index([1 3 4]) tmpoverlap]];
                end
            end
        end
    end
end

%create matriced with the data from the animals combined, instead of
%seperated in a cell matrix
ca1overlap2 = [];
ca3overlap2 = [];
for i = 1:length(ca1overlap)
    ca1overlap2 = [ca1overlap2; ca1overlap{i}];
    ca3overlap2 = [ca3overlap2; ca3overlap{i}];
end

%compare to output of pairwisescript
ca1 = indexmatch(ca1overlap2,ca1cmp,[1 2 3 4]);
ca3 = indexmatch(ca3overlap2,ca3cmp,[1 2 3 4]);
ca1data = [];
ca3data = [];
ca1data.index = ca1(:,1:4);
ca1data.overlap = ca1(:,5);
ca1data.rate1 = ca1(:,6);
ca1data.rate2 = ca1(:,7);
ca3data.index = ca3(:,1:4);
ca3data.overlap = ca3(:,5);
ca3data.rate1 = ca3(:,6);
ca3data.rate2 = ca3(:,7);

%plotting
ca1means = [];
ca3means = [];
ca1stderr = [];
ca3stderr = [];
for i = 1:length(ca1overlap)
    try
    ca1means(i) = mean(ca1overlap{i}(:,5));
    ca1stderr(i) = stderr(ca1overlap{i}(:,5));
    ca3means(i) = mean(ca3overlap{i}(:,5));
    ca3stderr(i) = stderr(ca3overlap{i}(:,5));
    end
end
figure
x = 1:6;
errorbar(x,ca1means(x),ca1stderr(x),'.');
hold on
errorbar(x,ca3means(x),ca3stderr(x),'r.');



ca1A = [ca1overlap{1}(:,5);ca1overlap{2}(:,5);ca1overlap{3}(:,5)];
ca1B = [ca1overlap{4}(:,5);ca1overlap{5}(:,5);ca1overlap{6}(:,5)];
ca3A = [ca3overlap{1}(:,5);ca3overlap{2}(:,5);ca3overlap{3}(:,5)];
ca3B = [ca3overlap{4}(:,5);ca3overlap{5}(:,5);ca3overlap{6}(:,5)];

