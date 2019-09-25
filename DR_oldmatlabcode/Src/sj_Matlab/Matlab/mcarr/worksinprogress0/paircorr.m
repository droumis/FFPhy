%Animal selection
animals = {'Bond','Frank','Ten'};

epochfilter = [];
epochfilter{1} = 'isequal($type,''run'')';

%Filter creation

cellpairfilter = {'allcomb', '(isequal($area, ''CA3'') && ($meanrate < 7))', '(isequal($area, ''CA3'') && ($meanrate < 7))'};
%cellpairfilter = {'allcomb', '$meanrate < 7', '$meanrate < 7'};

timefilter1 = {{'get2dstate', '((abs($velocity) > 1))'},{'getriptimes', '($nripples == 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'} }; %running, excluding ripples, eg place field activity
timefilter2 = {{'get2dstate', '((abs($velocity) < 4))'},{'getriptimes', '($nripples >= 1)', [], 'cellfilter', 'isequal($area, ''CA1'')'}}; %during ripples, bnot running

f1 = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetime', timefilter1);
f2 = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetime', timefilter2);

iterator = 'multicellanal2';

f1 = setfilteriterator(f1,iterator);
f2 = setfilteriterator(f2,iterator);

f1 = setfilterfunction(f1, 'calcpairxcorr', {'spikes','linpos'});
f2 = setfilterfunction(f2, 'calcpairxcorr', {'spikes','linpos'});

% run filter
f1 = runfilter(f1);
f2 = runfilter(f2);

%% Combine data and compare
paircount = 0;
pairs = [];
for animalnum = 1:length(f1) %for each animal
    tmpindex1 = zeros(length(f1(animalnum).output{1}),6); %length of number of pairs
    tmpindex2 = zeros(length(f2(animalnum).output{1}),6);
    for i = 1:length(f1(animalnum).output{1})
        tmpindex1(i,1:6) = f1(animalnum).output{1}(i).ind;
        tmpindex1(i,7) = i; %[day epoch tet cell tet cell pairnum]
    end
    for i = 1:length(f2(animalnum).output{1})
        tmpindex2(i,1:6) = f2(animalnum).output{1}(i).ind;
        tmpindex2(i,7) = i;
    end
    indexcolumns = [1 2 3 4 5 6];
    matchfind = indexmatch(tmpindex1,tmpindex2,indexcolumns); %[day epoch tet cell tet cell pair#index1 pair#index2]
    visited = zeros(size(matchfind,1),1);
    disp(animalnum)
    for j = 1:size(matchfind,1) %for each pair in both filter1 and filter2
        
        if (~visited(j))
            paircount = paircount+1;
            tmpmatchfind = matchfind(j,[1 3:6]);%[day tet cell tet cell]
            hits = ((matchfind(:,1)==tmpmatchfind(1)) & (matchfind(:,3)==tmpmatchfind(2)) & (matchfind(:,4)==tmpmatchfind(3)) & (matchfind(:,5)==tmpmatchfind(4)) & (matchfind(:,6)==tmpmatchfind(5)) );
                %find rows of matchfind that are included on the say day
                %but in different epochs
                
            %hits = ismember(matchfind(:,[1 3:6]),matchfind(j,[1 3:6]),'rows');
            pairs(paircount).index = matchfind(hits,1:6);
            corrind = matchfind(hits,7:8);
            for k = 1:size(corrind,1)
                %measure from f1
                
                pairs(paircount).filter1(pairs(paircount).index(k,2)).time = f1(animalnum).output{1}(corrind(k,1)).time;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).peakdist = f1(animalnum).output{1}(corrind(k,1)).peakdist;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).peakrate1 = f1(animalnum).output{1}(corrind(k,1)).peakrate1;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).peakrate2 = f1(animalnum).output{1}(corrind(k,1)).peakrate2;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).nspikes1 = f1(animalnum).output{1}(corrind(k,1)).nspikes1;
                pairs(paircount).filter1(pairs(paircount).index(k,2)).nspikes2 = f1(animalnum).output{1}(corrind(k,1)).nspikes2;
               
                %Add gamma phase
                pairs(paircount).filter1(pairs(paircount).index(k,2)).phase = f1(animalnum).output{1}(corrind(k,1)).phase;
                
                %measures from f2
                pairs(paircount).filter2(pairs(paircount).index(k,2)).time = f2(animalnum).output{1}(corrind(k,2)).time;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).nspikes1 = f2(animalnum).output{1}(corrind(k,2)).nspikes1;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).nspikes2 = f2(animalnum).output{1}(corrind(k,2)).nspikes2;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).peakdist = f2(animalnum).output{1}(corrind(k,1)).peakdist;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).peakrate1 = f2(animalnum).output{1}(corrind(k,1)).peakrate1;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).peakrate2 = f2(animalnum).output{1}(corrind(k,1)).peakrate2;
                
                %Add gamma phase
                pairs(paircount).filter2(pairs(paircount).index(k,2)).phase = f2(animalnum).output{1}(corrind(k,1)).phase;
               
            end              
            visited(hits) = 1;  %any unique pairs included on same day are set as "visited" so don't repeat.  this deals with same pairs in different epochs
        end
    end
    disp('done')
    
end

%% Compare second run ripple spike times with first run place fields
xlookup = -.5:.005:.5;
ylookup = 1:2:160;
glookup = -15*pi:pi/10:15*pi; 
ydensity = zeros(size(ylookup,2),1);
ydensity_phase = zeros(size(ylookup,2),1);
compare = [];
compare_phase = [];
for i = 1:length(pairs)
    try
        epoch = pairs(i).index(:,2);
        if length(epoch ) > 1
            if  ~isempty(pairs(i).filter1(epoch(1)).peakrate1) &&  ... 
                pairs(i).filter1(epoch(1)).peakrate1 > 5 && pairs(i).filter1(epoch(1)).peakrate2 > 5 && ...
                pairs(i).filter1(epoch(end)).nspikes1 < 50 && pairs(i).filter1(epoch(end)).nspikes2 < 50

                tmpspiketimes = pairs(i).filter2(epoch(2)).time';  %take cospike timing info from epoch 6
                peakdist = pairs(i).filter1(epoch(1)).peakdist; %get peakdist from epoch 2
                tmppeakdist = repmat(peakdist,length(tmpspiketimes),1);


                ydensity(lookup(pairs(i).filter1(epoch(1)).peakdist,ylookup)) = ydensity(lookup(peakdist,ylookup))+1;
                compare = [compare; [tmpspiketimes tmppeakdist]];

                tmpspiketimes = pairs(i).filter2(epoch(2)).phase';  %take cospike timing info from epoch 6
                peakdist = pairs(i).filter1(epoch(1)).peakdist; %get peakdist from epoch
                tmppeakdist = repmat(peakdist,length(tmpspiketimes),1);

                ydensity_phase(lookup(pairs(i).filter1(epoch(1)).peakdist,ylookup)) = ydensity_phase(lookup(peakdist,ylookup))+1;
                compare_phase = [compare_phase; [tmpspiketimes tmppeakdist]];
            end
        end
    end
end

sparseimage = [lookup(compare(:,2),ylookup) lookup(compare(:,1),xlookup)];
sparseimage(:,3) = 1;
sparseimage(:,3) = sparseimage(:,3)./(ydensity(sparseimage(:,1)));
x = spconvert(sparseimage);
x = full(x);
for i = 1:size(x,1)
    x(i,:) = x(i,:)/max(x(i,:));
end

figure
imagesc(xlookup,ylookup,x);
set(gca,'Ydir','normal','xtick',-0.5:0.5:0.5,'ytick',0:80:160)
ylabel('Distance between field peaks (cm)')
xlabel('Relative spike timing (s)')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_replay_time.png', m, d, y);
print('-dpng', savestring)

sparseimage = [lookup(compare_phase(:,2),ylookup) lookup(compare_phase(:,1),glookup)];
sparseimage(:,3) = 1;
sparseimage(:,3) = sparseimage(:,3)./(ydensity_phase(sparseimage(:,1)));
x_phase = spconvert(sparseimage);
x_phase = full(x_phase);
for i = 1:size(x_phase,1)
    x_phase(i,:) = x_phase(i,:)/max(x_phase(i,:));
end

figure
imagesc(glookup,ylookup,x_phase);
set(gca,'Ydir','normal','xtick',glookup(1):3*pi:glookup(end),'ytick',0:80:160)
ylabel('Distance between field peaks (cm)')
xlabel('Relative gamma phase of spikes (radians)')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_replay_gamma.png', m, d, y);
print('-dpng', savestring)

tmpcompare = compare(:,[2 1]); %compare: [spkiketimes, pkdist]
tmpcompare(:,3) = 1; %tmpcompare: [pkdist, spiketimes, ones]
tmpcompare = abs(tmpcompare); %mpcompare: [pkdist, spiketimes, ones]
[b,bint,r,rint,stats] = regress(tmpcompare(:,1),tmpcompare(:,2:3));

tmpcompare = compare_phase(:,[2 1]); %compare: [spkiketimes, pkdist]
tmpcompare(:,3) = 1; %tmpcompare: [pkdist, spiketimes, ones]
tmpcompare = abs(tmpcompare); %mpcompare: [pkdist, spiketimes, ones]
[bg,bintg,rg,rintg,statsg] = regress(tmpcompare(:,1),tmpcompare(:,2:3));

%R2 for time: 0.148, R2 for gamma:0.209
%Absolute value of residuals are significantly different, ttest, p<10e-10

%For CA3 only:
%R2 for time: 0.3352, R2 for gamma: 0.3547
%Absolute value of residuals are significantly different, ttest, p<10e-10

%For CA1 only:
%R2 for time: 0.2344, R2 for gamma: 0.312
%Absolute value of residuals are significantly different, ttest, p<10e-10