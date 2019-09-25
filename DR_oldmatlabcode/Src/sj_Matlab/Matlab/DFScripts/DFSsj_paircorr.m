

% Shantanu - Using the stock paircorr.m to plot reactivn peaks vs place
% field distance
% Shantanu - multicellanal is not working for me - So i am changing to
% singelcellanal like xcorrmeasures

% Match this to xcorrmeasures
% I will incorporate this into xcorrmeasures later


clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - 

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';
savefile = [savedir 'REGrp_paircorr'];
%savefile = [savedir 'REGrp_paircorr_traj'];



% If runscript, run Datafilter and save data
if runscript == 1
    
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'REc','REd','REe','REf'};
    %animals = {'REd','REe','REf'};
    %animals = {'REe'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    runepochfilter = 'isequal($type, ''run'')';
    sleepepochfilter = 'isequal($type, ''sleep'')';
    
    % Cell filter
    % -----------
    cellpairfilter = {'allcomb','strcmp($tag, ''PyrSR'')','strcmp($tag, ''PyrSR'')'};
    
    
    % Time filter
    % -----------
     
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % For peak distance during run - Needed only if linpos is used. Not for
    % linfields
    timefilter_place = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6},...
        {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };    
    
    % Note: Can use get2dstate for immobility time
    
    % Mainly For reactivn peaks during sleep. Using stimfilter as well   
%     timefiltersleep_onlyrip = {{'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3},...
%         {'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter}};
    timefiltersleep_onlyrip = {{'DFTFsj_getriptimes_nostim','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3}};
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    
    % No rip, run condition
    
    xrun = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
        cellpairfilter,'excludetime', timefilter_place,'iterator', iterator);
%     xrun = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
%         cellpairfilter,'iterator', iterator);
    
    % Rip, stop condition
    sleepcorr = createfilter('animal',animals,'days',dayfilter,'epochs',sleepepochfilter,'cellpairs',...
        cellpairfilter,'excludetime', timefiltersleep_onlyrip,'iterator', iterator);
    
    %f1 = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetime', timefilter1);
    %f2 = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetime', timefilter2);
    %f1 = setfilteriterator(f1,iterator);
    %f2 = setfilteriterator(f2,iterator);
    
     % Set analysis function
    % ----------------------
    xrun = setfilterfunction(xrun, 'DFAsj_peakdistance_linpos', {'spikes','linpos'});
    %xrun = setfilterfunction(xrun, 'DFAsj_peakdistance_traj', {'linfields'});
    sleepcorr = setfilterfunction(sleepcorr, 'DFAsj_calcpairxcorr', {'spikes'});
    
    disp('Done Filter Creation');
    
    % Run analysis
    % ------------
    xrun = runfilter(xrun);
    sleepcorr = runfilter(sleepcorr);
    
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript savedata
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

if ~exist('savedata')
    return
end




% -------------- Process Data --------------

useanim = 1:length(xrun); % Use all animals
usedays = [];   % Get for each animal separately

% Get run data - Only field peaks
for anidx = 1:length(useanim)
    an = useanim(anidx);
    for i = 1:length(xrun(an).output{1}), % This should be no. of values(pairs) X 2 epochs
        index{an}(i,:) = xrun(an).output{1}(i).index;
                
        allpeakdist{an}(i) = xrun(an).output{1}(i).peakdist;
        allpeakrate1{an}(i) = xrun(an).output{1}(i).peakrate1;
        allpeakrate2{an}(i) = xrun(an).output{1}(i).peakrate2;
    end
end

% Sleep corr - time of correlation peak
sleepcorrtime=[]; % Get timebase for sleepcorr only once
for anidx = 1:length(useanim)
    an = useanim(anidx);
    for i = 1:length(sleepcorr(an).output{1}), % This should be no. of values(pairs) X 3 epochs
        index_sl{an}(i,:) = sleepcorr(an).output{1}(i).index;
        
        corr_sl{an}(i,:) = sleepcorr(an).output{1}(i).c1vsc2;     
        normsmoothcorr_sl{an}(i,:) = sleepcorr(an).output{1}(i).normsmoothcorr; 
        nspikes1_sl{an}(i) = sleepcorr(an).output{1}(i).nspikes1;
        nspikes2_sl{an}(i) = sleepcorr(an).output{1}(i).nspikes2;
        time_sl{an}{i} = sleepcorr(an).output{1}(i).time;   
        maxtime_sl{an}{i} = sleepcorr(an).output{1}(i).maxtime;   
                
         if isempty(sleepcorrtime)
             if isfield(sleepcorr(an).output{1}(i),'timebase');
                 sleepcorrtime =  sleepcorr(an).output{1}(i).timebase;
             end
         end
    end
end


% Combine and calculate across animals

%xlookup = 0:.005:.5;
xlookup = -0.5:.005:.5;
ylookup = 0:4:200;
ydensity = zeros(size(ylookup,2),1);
compare = [];
maxcompare = [];

peakdistvec_ep4 = [];
peakdistvec_ep2 = [];
peakratevec1_ep4 = [];
peakratevec1_ep2 = [];
peakratevec2_ep4 = [];
peakratevec2_ep2 = [];


global_peakdist = [];
global_normsmoothcorr = [];


for anidx = 1:length(useanim)
    
    an = useanim(anidx);
    prefix = xrun(an).animal{1};
    animdirect = xrun(an).animal{2};
    if (animdirect(end) == '/')
        animdirect = animdirect(1:end-1);
    end
    
    % ******************
    % RUN DATA
    % Get idxs for two run epochs and parameters for curr animal. 
    curridxs_antot = index{an};
    ep2idxs = find(index{an}(:,2)==2);
    ep4idxs = find(index{an}(:,2)==4); % Ep2 and Ep4 MUST be same lengths
    curridxs_an = curridxs_antot(ep2idxs,:); % Going to collapse across run epochs   

    peakdist_ep2 = allpeakdist{an}(ep2idxs);
    peakdist_ep4 = allpeakdist{an}(ep4idxs);
    keep=find((peakdist_ep2./peakdist_ep4 >=0.5) & (peakdist_ep2./peakdist_ep4 <=2));
    peakdist_an = mean([peakdist_ep2(keep);peakdist_ep4(keep)]);
        
    peakrate1_ep2 = allpeakrate1{an}(ep2idxs(keep));  
    peakrate1_ep4 = allpeakrate1{an}(ep4idxs(keep));
    peakrate2_ep2 = allpeakrate2{an}(ep2idxs(keep));
    peakrate2_ep4 = allpeakrate2{an}(ep4idxs(keep));
    
    peakrate1_an = min([peakrate1_ep2; peakrate1_ep4]);
    peakrate2_an = min([peakrate2_ep2; peakrate2_ep4]);
    
    % Update indices
    ep2idxs = ep2idxs(keep);
    ep4idxs = ep4idxs(keep);
    curridxs_an = curridxs_an(keep,:);
    
    % ******************
    % SLEEP DATA
    % Get idxs for the three sleep epochs and params for current animal
    curridxs_antot_sl = index_sl{an};
    ep1idxs = find(index_sl{an}(:,2)==1);
    ep3idxs = find(index_sl{an}(:,2)==3);
    ep5idxs = find(index_sl{an}(:,2)==5); % Eps 1, 3 and 5 MUST be same lengths
    curridxs_an_sl = curridxs_antot_sl(ep1idxs,:); 
    
    % Update indices first
    ep1idxs = ep1idxs(keep);
    ep3idxs = ep3idxs(keep);
    ep5idxs = ep5idxs(keep);
    curridxs_an_sl = curridxs_an_sl(keep,:);
    
    nspikes1_ep5_an = nspikes1_sl{an}(ep5idxs);
    nspikes2_ep5_an = nspikes2_sl{an}(ep5idxs);
    
    % Get corrln time - put sleep epochs together?
    for i=1:length(ep5idxs)
        time1 = time_sl{an}{ep5idxs(i)};
        time2 = time_sl{an}{ep3idxs(i)};
        time_sl_ep5_an{i} = [time1 time2];
        %time_sl_ep5_an{i} = time_sl{an}{ep5idxs(i)};
    end
    
    
    % Get Nevents in raw correlation from -200ms to 200ms
    sw = 0.2;
    bins = find(abs(sleepcorrtime)<=sw);
        
    % Get entire normsmoothcorr
    for i=1:length(ep5idxs)
        
        rawcorr_ep5 = corr_sl{an}(ep5idxs(i),:);
        normsmoothcorr_ep5 = normsmoothcorr_sl{an}(ep5idxs(i),:);
        Neventscorr_ep5 = sum(rawcorr_ep5(bins));
        
        rawcorr_ep3 = corr_sl{an}(ep3idxs(i),:);
        normsmoothcorr_ep3 = normsmoothcorr_sl{an}(ep3idxs(i),:);
        Neventscorr_ep3 = sum(rawcorr_ep3(bins));
    
        if ((Neventscorr_ep3>=15) || (Neventscorr_ep5>=15))  % if either have >thresh spikes
            
            if ((Neventscorr_ep3<15) || (Neventscorr_ep5<15)) % if one has <thresh spikes, take the other
                if Neventscorr_ep3<15
                    normsmoothcorr_an(i,:) = normsmoothcorr_ep5;
                else
                    normsmoothcorr_an(i,:) = normsmoothcorr_ep3;
                end
            else % if both have >thresh spikes, take mean
                normsmoothcorr_an(i,:) = mean([normsmoothcorr_ep5; normsmoothcorr_ep3]);
            end
            
        else  % if neither have thresh spikes, take max
            
            if Neventscorr_ep3>Neventscorr_ep5
                normsmoothcorr_an(i,:) = normsmoothcorr_ep3;
            else
                normsmoothcorr_an(i,:) = normsmoothcorr_ep5;
            end          
            
        end
    
    end
    
        
    peakdistvec_ep4 = [peakdistvec_ep4,peakdist_ep4];
    peakdistvec_ep2 = [peakdistvec_ep2,peakdist_ep2];
    peakratevec1_ep4 = [peakratevec1_ep4,peakrate1_ep4];
    peakratevec1_ep2 = [peakratevec1_ep2,peakrate1_ep2];
    peakratevec2_ep4 = [peakratevec2_ep4,peakrate2_ep4];
    peakratevec2_ep2 = [peakratevec2_ep2,peakrate2_ep2];
    
    
    % Calculation: GATHER PARAMETERS
    
    for i=1:length(ep5idxs)  % All ep"n"idxs must be the same length
        
        if  ~isnan(peakrate1_an(i)) && ~isnan(peakrate2_an(i))  ...
                peakrate1_an(i) >= 3 && peakrate2_an(i) >= 3; %&& ...
                %nspikes1_ep5(i) < 50 && nspikes2_ep5(i) < 50
            
                tmpspiketimes = time_sl_ep5_an{i}';
                peakdist = peakdist_ep4(i); 
                tmppeakdist = repmat(peakdist,length(tmpspiketimes),1);
                
                ydensity(lookup(peakdist,ylookup)) = ydensity(lookup(peakdist,ylookup))+1;
                compare = [compare; [tmpspiketimes tmppeakdist]];
                
                % Max time - Time of Corrln Peak
                tmpspiketimes = maxtime_sl{an}{ep5idxs(i)}';
                tmppeakdist = repmat(peakdist,length(tmpspiketimes),1);                
                maxcompare = [maxcompare; [tmpspiketimes tmppeakdist]];
                
                % For taking actual normsmoothcorr
                global_peakdist = [global_peakdist; peakdist];
                global_normsmoothcorr = [global_normsmoothcorr; normsmoothcorr_an(i,:)];
            
        end
        
    end
   
end % end animals

       
% Based on corrln time
sparseimage = [lookup(compare(:,2),ylookup) lookup(compare(:,1),xlookup)];
sparseimage(find(sparseimage(:,2)==0),2) = 100; % Convert 0 if any to 1 - This shouldnt happen
sparseimage(:,3) = 1;
sparseimage(:,3) = sparseimage(:,3)./(ydensity(sparseimage(:,1)));
%sparseimage(:,3) = ydensity(sparseimage(:,1);
%sparseimage(find(sparseimage(:)<1))=1;
x = spconvert(sparseimage);
x = full(x);
for i = 1:size(x,1)
    x(i,:) = x(i,:)/max(x(i,:));
end

figure
imagesc(xlookup,ylookup,x);
set(gca,'Ydir','normal','xtick',-0.5:0.5:0.5,'ytick',0:40:160)
ylabel('Distance between field peaks (cm)')
xlabel('Relative spike timing (s)')
%colormap hot
 
% Based on normsmoothcorr
% yidxs = lookup(global_peakdist,ylookup);
% 
% for i=1:length(ylookup)
%     curr = find(yidxs==i);
%     if ~isempty(curr)
%         xn(i,:) = mean(global_normsmoothcorr(curr,:));
%     else
%         xn(i,:) = zeros(size(sleepcorrtime));
%     end
% end
%     
% for i = 1:size(xn,1)
%     xn(i,:) = xn(i,:)/max(xn(i,:));
% end
% figure
% imagesc(xlookup,ylookup,xn);
% set(gca,'Ydir','normal','xtick',-0.5:0.5:0.5,'ytick',0:40:160)
% ylabel('Distance between field peaks (cm)')
% xlabel('Relative spike timing (s)')

    



keyboard;

tmpcompare = compare(:,[2 1]); %compare: [spkiketimes, pkdist]
tmpcompare(:,3) = 1; %tmpcompare: [pkdist, spiketimes, ones]
tmpcompare = abs(tmpcompare); %mpcompare: [pkdist, spiketimes, ones]
[b,bint,r,rint,stats] = regress(tmpcompare(:,1),tmpcompare(:,2:3));



%% Plot median and inter-quartile range for gamma and time
ylookup = 0:10:160;
subs = lookup(compare(:,2),ylookup);

a = accumarray(subs,compare(:,1),[length(ylookup) 1],@(x) median(x));
a_seL = accumarray(subs,compare(:,1),[length(ylookup) 1],@(x) prctile(x,25));
a_seU = accumarray(subs,compare(:,1),[length(ylookup) 1],@(x) prctile(x,75));

figure
hold on
plot(ylookup,a,'r')
fill([ylookup ylookup(end:-1:1)],[a_seL; a_seU(end:-1:1)],'r','EdgeColor','none')
set(gca,'xtick',ylookup(1:2:end),'ytick',0:0.05:0.25,'ylim',[0 0.25],'xlim',[0 160])
xlabel('Distance between field peaks (cm)')
ylabel('Relative spike timing (sec)')

keyboard;
    


%% For which place field separation is time or gamma most predictive?
ylookup = [8 16 32 64];
time = nan(size(ylookup)); %phase = nan(size(ylookup));
subs = lookup(compare(:,2),ylookup);
for y = 1:length(ylookup)
    if sum(subs==y)>10
        time(y) = corr(compare(subs==y,1),compare(subs==y,2),'type','Spearman');
        %phase(y) = corr(compare_phase(subs==y,1),compare_phase(subs==y,2),'type','Spearman');
    end
end

% Run bootstrap analysis to compute confidence intervals on r
nboot = 1000;
qtime = nan(length(ylookup),nboot); %qphase = nan(length(ylookup),nboot);
qtime_all = nan(nboot,1); %qphase_all = nan(nboot,1);
t = cell(length(ylookup),1); %p = cell(length(ylookup),1);
for y = 1:length(ylookup)
    t{y} = compare(subs==y,:);
    %p{y} = compare_phase(subs==y,:);
end
for b = 1:nboot
    for y = 1:length(ylookup)
    	boot = t{y}(ceil(sum(subs==y)*rand(sum(subs==y),1)),:);
        qtime(y,b) = corr(boot(:,1),boot(:,2),'type','Spearman');    
    	%boot = p{y}(ceil(sum(subs==y)*rand(sum(subs==y),1)),:);
        %qphase(y,b) = corr(boot(:,1),boot(:,2),'type','Spearman');
    end
    boot = ceil(size(compare,1)*rand(size(compare,1),1));
    qtime_all(b) = corr(compare(boot,1),compare(boot,2),'type','Spearman');
    %qphase_all(b) = corr(compare_phase(boot,1),compare_phase(boot,2),'type','Spearman'); 
end
    
%qp = mean(qphase,2); se_p = std(qphase,[],2);
qt = mean(qtime,2); se_t = std(qtime,[],2);

figure
%bar(1:2:2*length(ylookup),qp,0.5,'r')
hold on
bar(2:2:2*length(ylookup),qt,0.5,'k')
%legend([{'Gamma phase'},{'Time'}],'Location','NorthWest')
legend([{'Time'}],'Location','NorthWest')
%errorbar2(1:2:2*length(ylookup),qp,se_p,'k')
errorbar2(2:2:2*length(ylookup),qt,se_t,'k')
set(gca,'xtick',1.5:2:2*length(ylookup),'xticklabel',ylookup,'xlim',[0 2*length(ylookup)+1])
xlabel('Distance between place field peaks')
ylabel('Spearman Correlation')

pvalue = zeros(length(ylookup),3);
%pvalue = zeros(length(ylookup),3);
for y = 1:length(ylookup)
    pvalue(y,1) = 1-sum(qtime(y,:)>0)./nboot;
    %pvalue(y,1) = 1-sum(qphase(y,:)>0)./nboot;
    %pvalue(y,2) = 1-sum(qtime(y,:)>0)./nboot;
    %pvalue(y,3) = 1-sum(qphase(y,:)>qt(y))./nboot;
end

% mean_all = [mean(qphase_all) mean(qtime_all)];
% se_all = [std(qphase_all) std(qphase_all)]; % Error; 2nd term must be qtime
% figure
% bar(1,mean_all(1),'r')
% hold on
% bar(2,mean_all(2),'k')
% errorbar2([1 2],mean_all,se_all,'k')
% set(gca,'xtick',[1 2],'xticklabel',[{'Gamma Phase'},{'Time'}],'xlim',[0.5 2.5])
% ylabel('Spearman correlation')
% 

mean_all = mean(qtime_all);
se_all = sem(qtime_all);
figure; hold on
bar(1,mean_all(1),'r')
%bar(2,mean_all(2),'k')
errorbar2([1],mean_all,se_all,'k')
set(gca,'xtick',[1],'xticklabel',{'Time'},'xlim',[0.2 1.8])
ylabel('Spearman correlation')










keyboard;





    
% ******************************************************
% ******************************************************
% ******************************************************
% ******************************************************
% ******************************************************


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
                %pairs(paircount).filter1(pairs(paircount).index(k,2)).phase = f1(animalnum).output{1}(corrind(k,1)).phase;
                
                %measures from f2
                pairs(paircount).filter2(pairs(paircount).index(k,2)).time = f2(animalnum).output{1}(corrind(k,2)).time;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).nspikes1 = f2(animalnum).output{1}(corrind(k,2)).nspikes1;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).nspikes2 = f2(animalnum).output{1}(corrind(k,2)).nspikes2;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).peakdist = f2(animalnum).output{1}(corrind(k,1)).peakdist;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).peakrate1 = f2(animalnum).output{1}(corrind(k,1)).peakrate1;
                pairs(paircount).filter2(pairs(paircount).index(k,2)).peakrate2 = f2(animalnum).output{1}(corrind(k,1)).peakrate2;
                
                %Add gamma phase
                %pairs(paircount).filter2(pairs(paircount).index(k,2)).phase = f2(animalnum).output{1}(corrind(k,1)).phase;
               
            end              
            visited(hits) = 1;  %any unique pairs included on same day are set as "visited" so don't repeat.  this deals with same pairs in different epochs
        end
    end
    disp('done')
    
end

%save('/data13/mcarr/RipplePaper/pairs.mat','pairs')









%% Compare sleep ripple spike times with preceding run place fields
xlookup = 0:.005:.5;
ylookup = 1:2:160;
%glookup = 0:pi/5:20*pi; 
ydensity = zeros(size(ylookup,2),1);
%ydensity_phase = zeros(size(ylookup,2),1);
compare = [];
%compare_phase = [];
for i = 1:length(pairs)
    try
        for epoch = [2 3; 4 5; 6 7]'
            if  ~isempty(pairs(i).filter1(epoch(1)).peakrate1) &&  ... 
                pairs(i).filter1(epoch(1)).peakrate1 > 3 && pairs(i).filter1(epoch(1)).peakrate2 > 3 && ...
                pairs(i).filter1(epoch(2)).nspikes1 < 50 && pairs(i).filter1(epoch(2)).nspikes2 < 50

                tmpspiketimes = pairs(i).filter2(epoch(2)).time';  %take cospike timing info from epoch n+1
                peakdist = pairs(i).filter1(epoch(1)).peakdist; %get peakdist from epoch n
                tmppeakdist = repmat(peakdist,length(tmpspiketimes),1);


                ydensity(lookup(pairs(i).filter1(epoch(1)).peakdist,ylookup)) = ydensity(lookup(peakdist,ylookup))+1;
                compare = [compare; [tmpspiketimes tmppeakdist]];

                tmpspiketimes = pairs(i).filter2(epoch(2)).phase';  %take cospike timing info from epoch 6
                peakdist = pairs(i).filter1(epoch(1)).peakdist; %get peakdist from epoch
                tmppeakdist = repmat(peakdist,length(tmpspiketimes),1);

                %ydensity_phase(lookup(pairs(i).filter1(epoch(1)).peakdist,ylookup)) = ydensity_phase(lookup(peakdist,ylookup))+1;
                %compare_phase = [compare_phase; [tmpspiketimes tmppeakdist]];
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
set(gca,'Ydir','normal','xtick',0:0.5:0.5,'ytick',0:80:160)
ylabel('Distance between field peaks (cm)')
xlabel('Relative spike timing (s)')
colormap hot


% Save figure
%[y, m, d] = datevec(date);
%savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_replay_time_sleep.png', m, d, y);
%print('-dpng', savestring)

%%%%%%% Phase %%%%%%%%
% sparseimage = [lookup(compare_phase(:,2),ylookup) lookup(compare_phase(:,1),glookup)];
% sparseimage(:,3) = 1;
% sparseimage(:,3) = sparseimage(:,3)./(ydensity_phase(sparseimage(:,1)));
% x_phase = spconvert(sparseimage);
% x_phase = full(x_phase);
% for i = 1:size(x_phase,1)
%     x_phase(i,:) = x_phase(i,:)/max(x_phase(i,:));
% end
% 
% figure
% imagesc(glookup,ylookup,x_phase);
% set(gca,'Ydir','normal','xtick',glookup(1):3*pi:glookup(end),'ytick',0:80:160)
% ylabel('Distance between field peaks (cm)')
% xlabel('Relative gamma phase of spikes (radians)')
% colormap hot
% 
% % Save figure
% [y, m, d] = datevec(date);
% savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_replay_gamma_sleep.png', m, d, y);
% print('-dpng', savestring)

tmpcompare = compare(:,[2 1]); %compare: [spkiketimes, pkdist]
tmpcompare(:,3) = 1; %tmpcompare: [pkdist, spiketimes, ones]
tmpcompare = abs(tmpcompare); %mpcompare: [pkdist, spiketimes, ones]
[b,bint,r,rint,stats] = regress(tmpcompare(:,1),tmpcompare(:,2:3));

%%%%%%% Phase %%%%%%%%
% tmpcompare = compare_phase(:,[2 1]); %compare: [spkiketimes, pkdist]
% tmpcompare(:,3) = 1; %tmpcompare: [pkdist, spiketimes, ones]
% tmpcompare = abs(tmpcompare); %mpcompare: [pkdist, spiketimes, ones]
% [bg,bintg,rg,rintg,statsg] = regress(tmpcompare(:,1),tmpcompare(:,2:3));

%R2 for time: 0.0386, R2 for gamma:0.0903
%Absolute value of residuals are significantly different, ranksum, p<10e-5
%For CA3 only:
%R2 for time: 0.111, R2 for gamma: 0.1105
%For CA1 only:
%R2 for time: 0.0595, R2 for gamma: 0.0656
%Absolute value of residuals are significantly different, ttest, p<10e-10


%% Plot median and inter-quartile range for gamma and time
ylookup = 0:10:160;
subs = lookup(compare(:,2),ylookup);

a = accumarray(subs,compare(:,1),[length(ylookup) 1],@(x) median(x));
a_seL = accumarray(subs,compare(:,1),[length(ylookup) 1],@(x) prctile(x,25));
a_seU = accumarray(subs,compare(:,1),[length(ylookup) 1],@(x) prctile(x,75));
%ap = accumarray(subs,compare_phase(:,1),[length(ylookup) 1],@(x) median(x)); 
%ap_seL = accumarray(subs,compare_phase(:,1),[length(ylookup) 1],@(x) prctile(x,25));
%ap_seU = accumarray(subs,compare_phase(:,1),[length(ylookup) 1],@(x) prctile(x,75));


figure
hold on
plot(ylookup,a,'r')
fill([ylookup ylookup(end:-1:1)],[a_seL; a_seU(end:-1:1)],'r','EdgeColor','none')
set(gca,'xtick',ylookup(1:2:end),'ytick',0:0.05:0.25,'ylim',[0 0.25],'xlim',[0 160])
xlabel('Distance between field peaks (cm)')
ylabel('Relative spike timing (sec)')

% Save figure
% [y, m, d] = datevec(date);
% savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_timevsdistance.pdf', m, d, y);
% print('-dpdf', savestring)

%%%%%%% Phase %%%%%%%%
% figure
% hold on
% plot(ylookup,ap,'r')
% fill([ylookup ylookup(end:-1:1)],[ap_seL; ap_seU(end:-1:1)],'r','EdgeColor','none')
% set(gca,'xtick',ylookup(1:2:end),'ytick',0:3*pi:15*pi,'ylim',[0 15*pi],'xlim',[0 160])
% xlabel('Distance between field peaks (cm)')
% ylabel('Relative gamma phase of spikes (radians)')
% % Save figure
% [y, m, d] = datevec(date);
% savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_gammavsdistance.pdf', m, d, y);
% print('-dpdf', savestring)


%% For which place field separation is time or gamma most predictive?
ylookup = [8 16 32 64];
time = nan(size(ylookup)); %phase = nan(size(ylookup));
subs = lookup(compare(:,2),ylookup);
for y = 1:length(ylookup)
    if sum(subs==y)>10
        time(y) = corr(compare(subs==y,1),compare(subs==y,2),'type','Spearman');
        %phase(y) = corr(compare_phase(subs==y,1),compare_phase(subs==y,2),'type','Spearman');
    end
end

% Run bootstrap analysis to compute confidence intervals on r
nboot = 1000;
qtime = nan(length(ylookup),nboot); %qphase = nan(length(ylookup),nboot);
qtime_all = nan(nboot,1); %qphase_all = nan(nboot,1);
t = cell(length(ylookup),1); %p = cell(length(ylookup),1);
for y = 1:length(ylookup)
    t{y} = compare(subs==y,:);
    %p{y} = compare_phase(subs==y,:);
end
for b = 1:nboot
    for y = 1:length(ylookup)
    	boot = t{y}(ceil(sum(subs==y)*rand(sum(subs==y),1)),:);
        qtime(y,b) = corr(boot(:,1),boot(:,2),'type','Spearman');
    
    	%boot = p{y}(ceil(sum(subs==y)*rand(sum(subs==y),1)),:);
        %qphase(y,b) = corr(boot(:,1),boot(:,2),'type','Spearman');
    end
    boot = ceil(size(compare,1)*rand(size(compare,1),1));
    qtime_all(b) = corr(compare(boot,1),compare(boot,2),'type','Spearman');
    %qphase_all(b) = corr(compare_phase(boot,1),compare_phase(boot,2),'type','Spearman'); 
end
    
%qp = mean(qphase,2); se_p = std(qphase,[],2);
qt = mean(qtime,2); se_t = std(qtime,[],2);

figure
%bar(1:2:2*length(ylookup),qp,0.5,'r')
hold on
bar(2:2:2*length(ylookup),qt,0.5,'k')
%legend([{'Gamma phase'},{'Time'}],'Location','NorthWest')
legend([{'Time'}],'Location','NorthWest')
%errorbar2(1:2:2*length(ylookup),qp,se_p,'k')
errorbar2(2:2:2*length(ylookup),qt,se_t,'k')
set(gca,'xtick',1.5:2:2*length(ylookup),'xticklabel',ylookup,'xlim',[0 2*length(ylookup)+1])
xlabel('Distance between place field peaks')
ylabel('Spearman Correlation')

pvalue = zeros(length(ylookup),3);
%pvalue = zeros(length(ylookup),3);
for y = 1:length(ylookup)
    pvalue(y,1) = 1-sum(qtime(y,:)>0)./nboot;
    %pvalue(y,1) = 1-sum(qphase(y,:)>0)./nboot;
    %pvalue(y,2) = 1-sum(qtime(y,:)>0)./nboot;
    %pvalue(y,3) = 1-sum(qphase(y,:)>qt(y))./nboot;
end

% mean_all = [mean(qphase_all) mean(qtime_all)];
% se_all = [std(qphase_all) std(qphase_all)]; % Error; 2nd term must be qtime
% figure
% bar(1,mean_all(1),'r')
% hold on
% bar(2,mean_all(2),'k')
% errorbar2([1 2],mean_all,se_all,'k')
% set(gca,'xtick',[1 2],'xticklabel',[{'Gamma Phase'},{'Time'}],'xlim',[0.5 2.5])
% ylabel('Spearman correlation')
% 

mean_all = mean(qtime_all);
se_all = sem(qtime_all);
figure
bar(1,mean_all(1),'r')
hold on
%bar(2,mean_all(2),'k')
errorbar2([1],mean_all,se_all,'k')
set(gca,'xtick',[1],'xticklabel',{'Time'},'xlim',[0.2 1.8])
ylabel('Spearman correlation')






% % Save figure
% [y, m, d] = datevec(date);
% savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_gammavstime_correlation_sleep.pdf', m, d, y);
% print('-dpdf', savestring)

%FOR RUN: gamma is significantly correlated with distance between
%placefield peaks for last three groups, p<0.01 p<0.001 p<0.001

%time is significantly correlated with distance between placefield peaks
%for last two groups, p<0.001

%gamma is significantly more correlated with time for place fields with
%distances far apart, p<0.001

%FOR SLEEP: Gamma and time are significantly correlated with distance
%between place field peaks for all 4 groups, p<1e-5, gamma is significantly
%more correlated than time for first and last groups, p<1e-5





% %% COMPARE SLEEP AND RUN PAIR CORRELATIONS
% run_compare = compare; run_phase = compare_phase;
% sleep_compare = compare; sleep_phase = compare_phase;
% 
% q = [corr(run_compare(:,1),run_compare(:,2),'type','Spearman')...
%     corr(sleep_compare(:,1),sleep_compare(:,2),'type','Spearman')...
%     corr(run_phase(:,1),run_phase(:,2),'type','Spearman')...
%     corr(sleep_phase(:,1),sleep_phase(:,2),'type','Spearman')];
% 
% % Run bootstrap analysis to compute confidence intervals on r
% nboot = 1000;
% qboot = nan(nboot,4);
% for b = 1:nboot
%   	boot = ceil(size(run_compare,1)*rand(size(run_compare,1),1));
%     qboot(b,1) = corr(run_compare(boot,1),run_compare(boot,2),'type','Spearman');
%     qboot(b,3) = corr(run_phase(boot,1),run_phase(boot,2),'type','Spearman');
%   
%     boot = ceil(size(sleep_compare,1)*rand(size(sleep_compare,1),1));
%     qboot(b,2) = corr(sleep_compare(boot,1),sleep_compare(boot,2),'type','Spearman');
%     qboot(b,4) = corr(sleep_phase(boot,1),sleep_phase(boot,2),'type','Spearman');
% end
%     
% 
% figure
% bar(1:4,mean(qboot(:,[1 3 2 4])),1,'r')
% hold on
% errorbar2(1:4,mean(qboot(:,[1 3 2 4])),std(qboot(:,[1 3 2 4])),'k')
% set(gca,'xtick',1:1:4,'xticklabel',[{'Time'},{'Phase'},{'Time'},{'Phase'}])
% ylabel('Spearman Correlation')
% 
% pvalue = zeros(6,1);
% pvalue(1:4) = 1-sum(qboot>0)./nboot;
% pvalue(5) = sum(qboot(1,:)>q(2))./nboot;
% pvalue(6) = sum(qboot(3,:)>q(4))./nboot;
% 
% % % Save figure
% [y, m, d] = datevec(date);
% savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_gammavstime_runvssleep.pdf', m, d, y);
% print('-dpdf', savestring)
% 
% %Behavioral sessions are significantly greater than rest sessions, p<0.001
% %Relative gamma phase is greater than relative spike timing p<0.01















%%%%% OTHER CALCS AND PLOTS %%%%%%%%%

% %% Compare remote run ripple spike times with run place fields
% xlookup = 0:.005:.5;
% ylookup = 1:2:160;
% glookup = 0:pi/5:20*pi; 
% ydensity = zeros(size(ylookup,2),1);
% ydensity_phase = zeros(size(ylookup,2),1);
% compare = [];
% compare_phase = [];
% for i = 1:length(pairs)
%     try
%         for epoch = [2 4; 2 6; 4 6]'
%             if  ~isempty(pairs(i).filter1(epoch(1)).peakrate1) &&  ... 
%                 pairs(i).filter1(epoch(1)).peakrate1 > 3 && pairs(i).filter1(epoch(1)).peakrate2 > 3 && ...
%                 pairs(i).filter1(epoch(2)).nspikes1 < 50 && pairs(i).filter1(epoch(2)).nspikes2 < 50
%                 
%                 tmpspiketimes = pairs(i).filter2(epoch(2)).time';  %take cospike timing info from epoch 6
%                 peakdist = pairs(i).filter1(epoch(1)).peakdist; %get peakdist from epoch 2
%                 tmppeakdist = repmat(peakdist,length(tmpspiketimes),1);
% 
% 
%                 ydensity(lookup(pairs(i).filter1(epoch(1)).peakdist,ylookup)) = ydensity(lookup(peakdist,ylookup))+1;
%                 compare = [compare; [tmpspiketimes tmppeakdist]];
% 
%                 tmpspiketimes = pairs(i).filter2(epoch(2)).phase';  %take cospike timing info from epoch 6
%                 peakdist = pairs(i).filter1(epoch(1)).peakdist; %get peakdist from epoch
%                 tmppeakdist = repmat(peakdist,length(tmpspiketimes),1);
% 
%                 ydensity_phase(lookup(pairs(i).filter1(epoch(1)).peakdist,ylookup)) = ydensity_phase(lookup(peakdist,ylookup))+1;
%                 compare_phase = [compare_phase; [tmpspiketimes tmppeakdist]];
%                
%             end
%         end
%     end
% end
% compare = abs(compare);
% sparseimage = [lookup(compare(:,2),ylookup) lookup(compare(:,1),xlookup)];
% sparseimage(:,3) = 1;
% sparseimage(:,3) = sparseimage(:,3)./(ydensity(sparseimage(:,1)));
% x = spconvert(sparseimage);
% x = full(x);
% for i = 1:size(x,1)
%     x(i,:) = x(i,:)/max(x(i,:));
% end
% 
% figure
% imagesc(xlookup,ylookup,x);
% set(gca,'Ydir','normal','xtick',xlookup(1):0.1:xlookup(end),'ytick',0:20:160)
% ylabel('Distance between field peaks (cm)')
% xlabel('Relative spike timing (s)')
% colormap hot
% 
% % % Save figure
% [y, m, d] = datevec(date);
% savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_replay_time_run.png', m, d, y);
% print('-dpng', savestring)
% 
% compare_phase = abs(compare_phase);
% sparseimage = [lookup(compare_phase(:,2),ylookup) lookup(compare_phase(:,1),glookup)];
% sparseimage(:,3) = 1;
% sparseimage(:,3) = sparseimage(:,3)./(ydensity_phase(sparseimage(:,1)));
% x_phase = spconvert(sparseimage);
% x_phase = full(x_phase);
% for i = 1:size(x_phase,1)
%     x_phase(i,:) = x_phase(i,:)/max(x_phase(i,:));
% end
% 
% figure
% imagesc(glookup,ylookup,x_phase);
% set(gca,'Ydir','normal','xtick',glookup(1):3*pi:glookup(end),'ytick',0:20:160)
% ylabel('Distance between field peaks (cm)')
% xlabel('Relative gamma phase of spikes (radians)')
% colormap hot
% % 
% % % Save figure
% [y, m, d] = datevec(date);
% savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_pairwise_replay_gamma_run.png', m, d, y);
% print('-dpng', savestring)
% 
% tmpcompare = compare(:,[2 1]); %compare: [spkiketimes, pkdist]
% tmpcompare(:,3) = 1; %tmpcompare: [pkdist, spiketimes, ones]
% %tmpcompare = abs(tmpcompare); %mpcompare: [pkdist, spiketimes, ones]
% [b,bint,r,rint,stats] = regress(tmpcompare(:,1),tmpcompare(:,2:3));
% 
% tmpcompare = compare_phase(:,[2 1]); %compare: [spkiketimes, pkdist]
% tmpcompare(:,3) = 1; %tmpcompare: [pkdist, spiketimes, ones]
% %tmpcompare = abs(tmpcompare); %mpcompare: [pkdist, spiketimes, ones]
% [bg,bintg,rg,rintg,statsg] = regress(tmpcompare(:,1),tmpcompare(:,2:3));
% 
% %R2 for time: 0.1405, R2 for gamma:0.2117
% %Absolute value of residuals are significantly different, ranksum, p<0.001
% 
% %For CA3 only:
% %R2 for time: 0.2537, R2 for gamma: 0.3196
% %Absolute value of residuals are significantly different, ttest, p<10e-10
% 
% %For CA1 only:
% %R2 for time: 0.123, R2 for gamma: 0.1811
% %Absolute value of residuals are significantly different, ttest, p<10e-10










