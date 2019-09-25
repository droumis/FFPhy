%Animal selection
%-----------------------------------------------------
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};

%animals = {'Miles','Nine','Ten'};
animals = {'Conley','Bond','Frank'};
%animals = {'Conley'};
%animals = {'Bond'};
%animals = {'Frank'};

%animals = {'Miles','Conley','Ten','Bond','Frank','Nine'};
%animals = {'Dudley','Miles','Conley','Alex','Bond','Frank','Nine'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------

epochfilter = [];
epochfilter{1} = ['isequal($type,''sleep'')|isequal($type,''run'')'];
%epochfilter{1} = ['($dailyexposure == 1) & isequal($environment, ''TrackB'')'];  


f1 = [];
f2 = [];
cmp = [];



%Filter creation
%--------------------------------------------------------

dayfilter = ''; %leave blank to take all days

%cellpairfilter = {'allcomb', '(isequal($area, ''CA1'') && ($meanrate < 7))', '(isequal($area, ''CA1'') && ($meanrate < 7))'};
cellpairfilter = {'allcomb', '(($meanrate < 7))', '(($meanrate < 7))'};
%cellpairfilter = {'difftet', '(($meanrate < 7))', '(($meanrate < 7))'};


%timefilter1 = {{'get2dstate', '((abs($velocity) > 1))'},{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',2} };
%timefilter2 = {{'get2dstate', '((abs($velocity) < 1))'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};

timefilter1 = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',2} };
%timefilter2 = {{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};

timefilter2 = {{'get2dstate', '($immobilitytime > 0)&($immobilitytime < 5)'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};
%timefilter2 = {{'getsplittimes','$includeseg == 1',2,2},{'get2dstate', '($immobilitytime > 0)&($immobilitytime < 5)'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};


%timefilter2 = {{'gethighthetatimes', '($nhightheta == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', '(isequal($area, ''CA1''))|(isequal($area, ''CA3''))','minstd',3}};
%timefilter2 = {};


%timefilter2 = {{'get2dstate', '((abs($velocity) < 3))'}};
%timefilter2 = {{'getlinstate', '($distance2well < 1.5)',6},{'getriptimes', '($nripples > 2)', [], 'cellfilter', '(isequal($area, ''CA1''))','minstd',3}};

f1 = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetime', timefilter1);
f2 = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetime', timefilter2);

%-----------------------------------------------------------

%Pairwise comparison
%-----------------------------------------------------------
%iterator = 'singlecellanal';
iterator = 'multicellanal2';


f1 = setfilteriterator(f1,iterator);
f2 = setfilteriterator(f2,iterator);


f1 = setfilterfunction(f1, 'calcpairxcorr', {'spikes','linpos'});
f2 = setfilterfunction(f2, 'calcpairxcorr', {'spikes','linpos'});



%f1 = setfilterfunction(f1, 'calcrippleactivationprob', {'spikes','ripples','cellinfo'},'appendindex',1,'ratio',1);
%f2 = setfilterfunction(f2, 'calcrippleactivationprob', {'spikes','ripples','cellinfo'},'appendindex',1,'ratio',1);

f1 = runfilter(f1);
f2 = runfilter(f2);
paircount = 0;
pairs = [];
for animalnum = 1:length(f1)
    tmpindex1 = zeros(length(f1(animalnum).output{1}),6);
    tmpindex2 = zeros(length(f2(animalnum).output{1}),6);
    for i = 1:length(f1(animalnum).output{1})
        tmpindex1(i,1:6) = f1(animalnum).output{1}(i).ind;
        tmpindex1(i,7) = i;
    end
    for i = 1:length(f2(animalnum).output{1})
        tmpindex2(i,1:6) = f2(animalnum).output{1}(i).ind;
        tmpindex2(i,7) = i;
    end
    indexcolumns = [1 2 3 4 5 6];
    matchfind = indexmatch(tmpindex1,tmpindex2,indexcolumns);
    visited = zeros(size(matchfind,1),1);
    disp(animalnum)
    for j = 1:size(matchfind,1)
        
        if (~visited(j))
            paircount = paircount+1;
            tmpmatchfind = matchfind(j,[1 3:6]);
            hits = ((matchfind(:,1)==tmpmatchfind(1)) & (matchfind(:,3)==tmpmatchfind(2)) & (matchfind(:,4)==tmpmatchfind(3)) & (matchfind(:,5)==tmpmatchfind(4)) & (matchfind(:,6)==tmpmatchfind(5)) );
            
            %hits = ismember(matchfind(:,[1 3:6]),matchfind(j,[1 3:6]),'rows');
            pairs(paircount).index = matchfind(hits,1:6);
            corrind = matchfind(hits,7:8);
            for k = 1:size(corrind,1)
                pairs(paircount).filter1(pairs(paircount).index(k,2)).time = f1(animalnum).output{1}(corrind(k,1)).time;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).peakdist = f1(animalnum).output{1}(corrind(k,1)).peakdist;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).peakrate1 = f1(animalnum).output{1}(corrind(k,1)).peakrate1;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).peakrate2 = f1(animalnum).output{1}(corrind(k,1)).peakrate2;
                %pairs(paircount).filter1(pairs(paircount).index(k,2)).crosscorr = f1(animalnum).output{1}(corrind(k,1)).c1vsc2;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).nspikes1 = f1(animalnum).output{1}(corrind(k,1)).nspikes1;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).nspikes2 = f1(animalnum).output{1}(corrind(k,1)).nspikes2;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).time = f2(animalnum).output{1}(corrind(k,2)).time;
                %pairs(paircount).filter2(pairs(paircount).index(k,2)).crosscorr = f2(animalnum).output{1}(corrind(k,2)).c1vsc2;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).nspikes1 = f2(animalnum).output{1}(corrind(k,2)).nspikes1;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).nspikes2 = f2(animalnum).output{1}(corrind(k,2)).nspikes2;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).peakdist = f2(animalnum).output{1}(corrind(k,1)).peakdist;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).peakrate1 = f2(animalnum).output{1}(corrind(k,1)).peakrate1;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).peakrate2 = f2(animalnum).output{1}(corrind(k,1)).peakrate2;
            end              
            visited(hits) = 1;
        end
    end
    disp('done')
    
end

xlookup = [-.5:.005:.5];
%xlookup = [0:.005:.5];
ylookup = [1:2:300];
%ylookup = [1:20:300];
ydensity = zeros(size(ylookup,2),1);
compare = [];
for i = 1:length(pairs)
    try
        if ( (~isempty(pairs(i).filter1(6).peakrate1)) &&  ... 
                (pairs(i).filter1(2).peakrate1 > 1) && (pairs(i).filter1(2).peakrate2 > 1) && ...            
             ((pairs(i).filter1(3).nspikes1 < 50) && (pairs(i).filter1(3).nspikes2 < 50)) && ...
             (pairs(i).index(1,3) ~= pairs(i).index(1,5)))
%            totalspikes = pairs(i).filter1(2).nspikes1*pairs(i).filter1(2).nspikes2;
%            tmppeakdist = pairs(i).filter1(2).peakdist;
%            cbins = find(abs(pairs(i).filter1(2).time < .05));
%            cofirespikes = sum(pairs(i).filter1(2).crosscorr(cbins));
%            [mval,peakind] = max(pairs(i).filter2(3).crosscorr+pairs(i).filter2(5).crosscorr);
%            %[mval,peakind] = max(pairs(i).filter2(3).crosscorr);
%            peakoffset = abs(pairs(i).filter2(3).time(peakind));
%            %cbins2 = find(abs(pairs(i).filter2(6).time < .05));
%            %cofirespikes2 = sum(pairs(i).filter2(6).crosscorr(cbins2));
%            %totalspikes2 = pairs(i).filter2(6).nspikes1*pairs(i).filter2(6).nspikes2;
%            
%            if (mval > 10)
%             
%             compare = [compare;[tmppeakdist peakoffset]];
%             %compare = [compare; [(cofirespikes/totalspikes) (cofirespikes2/totalspikes2) mval]];
%             
%            end
           
            tmpspiketimes = abs([pairs(i).filter2(3).time'; pairs(i).filter2(5).time']); 
            
            %tmpspiketimes = abs([pairs(i).filter2(6).time']);
            peakdist = pairs(i).filter1(2).peakdist;
            tmppeakdist = repmat(peakdist,length(tmpspiketimes),1);
            
            
            ydensity(lookup(pairs(i).filter1(2).peakdist,ylookup)) = ydensity(lookup(peakdist,ylookup))+1;
            compare = [compare; [tmpspiketimes tmppeakdist]];
           
        end
    
    end
    
end






%tmphist = zeros(length(ylookup),length(xlookup));
sparseimage = [lookup(compare(:,2),ylookup) lookup(compare(:,1),xlookup)];
%ydensity = rowcount(ylookup',sparseimage(:,1));
sparseimage(:,3) = 1;
sparseimage(:,3) = sparseimage(:,3)./(ydensity(sparseimage(:,1)));
x = spconvert(sparseimage);
x = full(x);
for i = 1:size(x,1)
    x(i,:) = x(i,:)/max(x(i,:));
end

%x(find(x>2)) = 2;
% g = gaussian2(1,(10*1));
% x = filter2(g,x);

%[x, bx, by] = hist2(compare,[-.5:.02:.5],[0:2:250]);
figure
imagesc(x);
set(gca,'Ydir','normal');





tmpcompare = compare(:,[2 1]);
tmpcompare(:,3) = 1;
[b,bint,r,rint,stats] = regress(tmpcompare(:,1),tmpcompare(:,2:3));






for i = 1:length(pairs)
    try
        if ((pairs(i).filter1(2).nspikes1 > 100) && (pairs(i).filter1(2).nspikes2 > 100))
            subplot(3,1,1)
            bar(pairs(i).filter2(1).time,pairs(i).filter2(1).crosscorr);
            subplot(3,1,2)
            bar(pairs(i).filter1(2).time,(pairs(i).filter1(2).crosscorr+pairs(i).filter1(4).crosscorr));
            subplot(3,1,3)
            bar(pairs(i).filter2(3).time,(pairs(i).filter2(3).crosscorr+pairs(i).filter2(5).crosscorr));
            pause
        end
    catch
        disp('error')
    end
end

