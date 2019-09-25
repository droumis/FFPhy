%% Here, I look at two main things: 1) Does the presence of an auditory stimulus affect ripples? 2) What is going on, looking at EEG data, in the three brain regions during stimuli? 

% I have presented these results in my rotation ending presentation.  See
% file: "Rotation Final Presentation Winter 2013_3.ppt" The substance has not been
% changed since then. 

% I noted results here in powerpoint: "Borg - 1st analysis.odp"

% Everything here is for Borg, day 17, epoch 1.  It is not very complete or
% clean.  Much statistical analysis is needed.


%% PART I - STIMULUS EFFECT ON RIPPLES.  First, extract ripple events

% EXTRACTING SWR TIMES
% sj_rippledayprocess('/opt/data15/gideon/Brg','brg',17) %Do this to
% filter the raw EEG with ripple band.  

%EXTRACT RIPPLE EVENTS
extractripplesnopos('/home/lucas/gideon_data_copy/data15/gideon/Brg', 'brg', 17, -1, 0.015, 5) %the last argument is the STD threshold

%Loading ripples from all tetrodes.
allRipInd = cell(1,21);

load('/home/lucas/gideon_data_copy/data15/gideon/Brg/brgripples17.mat')

for tet=1:21,
    allRipInd{tet}=ripples{17}{1}{tet}.startind; % The cell dimensions: {day}{epoch}{tetrode}
end
hold on; 

ripInds=allRipInd{10}; %Will use tet 10 for SWR events as triggers, as LFP and event timepoints looked good.


%% Second, extram auditory stimulus times.  

%the output allTrialsPerStimInd gives the index of each stim, organized as
%a matrix of 75 (stim) X 20 (trials)
% There are 15 frequencies and 5 attenuations, so that the order of stim is 20 repeats of (freq1-att1), then 20 repeats of (freq1-att2), etc
% Organized in 20x75 array as: rows: trial# increases going down.  columns going right: freq#1-atten#1 f1-a2 ... f2-a1 ...
% you can look in getResponsesLFP to see how to extract
[allTrialsPerStim allTrialsPerStimInd stimEEGtimes stimInds]=findStimBrg('/opt/data15/gideon/Brg/EEG/Brgeeg', '17','1','22','/home/gideon/Code/Matlab/makingSounds/soundSeries/toneSweep1/freqsAttensToneSweep1.mat',0);


%% Third: Does stimulus affect ripple rates?  Find the stim-triggered average of SWR events.  Additionally, create a peri-stim histograms of ripples relative to stimulus.

% finding the stim-triggered number of SWR events in time bins in each stim's window.   
stimTrigRip_post = NaN(length(stimInds),2);
stimTrigRipBefore = NaN(length(stimInds),2);

%enter parameters
timeWindowInS=3;
for m=1:length(stimSorted_vector)
samplesPerWindow=1500*timeWindowInS;

%enter histogram specific parameters
bin_sec= 0.15;  %bin size in sec
bin=bin_sec*1500;
periStimRip_all=NaN(1+samplesPerWindow/bin,length(stimInds)); %matrix to hold all stims' windows.


%calculate histogram, and the average values of # of ripples before and
%after stim. 
for n=1:length(stimInds);
stimTrig=stimInds(n);

stimTrigWindowEnd=stimTrig+samplesPerWindow/2;
stimTrigWindowStart=stimTrig-samplesPerWindow/2;

howManyRipplesAfterStim=length(find(ripInds>stimTrig & ripInds<stimTrigWindowEnd)); % finds the number of ripple events that had indices between low and high defined for stim window, puts in array.
stimTrigRip_post(n,:)=[stimTrig howManyRipplesAfterStim]; % Gives array with Stim Indices and corresponding # of SWR in post-stim window

numRipBeforeStim=length(find(ripInds>stimTrigWindowStart & ripInds<stimTrig)); %same, but for rips before stim. 
stimTrigRipBefore(n,:)=[stimTrig numRipBeforeStim];

stim_edges=stimTrigWindowStart:bin:stimTrigWindowEnd; %edges of the histogram bins.
periStimRip_all(:,n)=histc(ripInds,stim_edges); %calculating array of histograms for each stim (stim=column)
end

%plotting histogram averaged over all Stim windows
periStimRip_mean=mean(periStimRip_all, 2); %use this to plot the mean per bin
periStimRip_sum=sum(periStimRip_all,2); % or have the histogram show total # of rips per bin.

figure
bar(0:bin_sec:timeWindowInS,periStimRip_sum,'histc'), axis tight, hold on, plot(timeWindowInS/2,0,'r'); title('Peri-stimulus ripple events (all stimuli in an epoch)')

% average values for peri-stim rips - show these raw values.  
stimTrigRip_mean=mean(stimTrigRip_post(:,2)) %post stim
stimTrigRip_std=std(stimTrigRip_post(:,2))
stimTrigRipBefore_mean=mean(stimTrigRipBefore(:,2)) %pre stim
stimTrigRipBefore_std=std(stimTrigRipBefore(:,2))

%% Fourth, plotting histograms for specific auditory stimuli (i.e. sort stimuli by auditory frequencies, attenuations, or trials)

%histogram inputs - which stim do you want to see?  Attn # is loudest at 5.
trials=1:20;
freqs=1:15;
attns=3:5;

%array with sampling indices for values.  See all the way up in this doc
%for details.  In short, each Frequency/attenuation combo was give one
%index number in the same dimension (trial# takes up the second dimension).  The better was would have been to have them take up
%separate dimensions.  

FrAtInd_list=[]; %The index numbers for your desired freq/attn combos.
clear n; clear m;
for n=1:length(freqs);
    fr=freqs(n);
    for at=attns
        FrAtInd=(fr-1)*5+at;
        FrAtInd_list=[FrAtInd_list FrAtInd];
    end
end

sort(FrAtInd_list)

% Row indices for freq-attn sorted stimuli. 
stimSorted=allTrialsPerStimInd(trials,FrAtInd_list);
stimSorted_vector=reshape(stimSorted,1,numel(stimSorted));

periStimRip_sorted=NaN(1+samplesPerWindow/bin,length(stimSorted_vector));

% takes all the stimuli you want and spits out an array containing a
% histogram (telling you how many ripples in each bin) over all windows of those stimuli
for m=1:length(stimSorted_vector)
    stimTrigSorted=stimSorted_vector(m);
    stimTrigWindowStartSorted=stimTrigSorted-samplesPerWindow/2;
    stimTrigWindowEndSorted=stimTrigSorted+samplesPerWindow/2;
    stim_edgesSorted=stimTrigWindowStartSorted:bin:stimTrigWindowEndSorted;
    periStimRip_sorted(:,m)=histc(ripInds,stim_edgesSorted);  
end


%% Fifth, AVERAGE # OF RIPPLES IN STIM OR NON-STIM PERIODS
% Here, I calculate single values for the mean # of ripples in various
% stim, no stim epochs.  Use these to compare to the histogram.

%rate of ripples in stim period
TotalNoRipsStim=length(find(ripInds<max(stimInds) & ripInds>min(stimInds))) %total No of ripples during stim period (i.e. period actually containing stimuli)
TotalStimDur=(max(stimInds)-min(stimInds))/1500 %total stim duration in seconds
RipRateDurStim_mean=TotalNoRipsStim/TotalStimDur %rate of ripples during stim period, rips/sec


% rate of ripples in the non-Stimulated epoch (i.e. after stim presentation ended)
noStimEpochMin= max(stimInds)+10*1500; %start 10 sec after last sound
noStimEpochMax=length(allEEGs(:,10))-5*1500; % end of EEG data, minus 5 sec.    
numNoStimRip=length(find(ripInds>noStimEpochMin & ripInds<noStimEpochMax)); %gives number of ripples over the non-stim epoch
noStimEpoch_meanRip=numNoStimRip/((noStimEpochMax-noStimEpochMin)/1500) %gives average Rip for the entire epoch. per sec., 

%Average # of stimuli per bin size (in both stim and no stim periods), to plot with peri-stim histogram as control

StimPeriodRipsInBin = RipRateDurStim_mean*bin_sec
NoStimPeriodRipsInBin = noStimEpoch_meanRip*bin_sec

%plotting the histogram again, but this time indicating the mean rate, for
%comparison.
periStimRipSorted_mean=mean(periStimRip_sorted,2);
periStimRipSorted_sum=sum(periStimRip_sorted,2);
mean(periStimRipSorted_mean)

figure
bar(0:bin_sec:timeWindowInS,periStimRipSorted_sum,'histc'), hold on, plot(timeWindowInS/2,0,'r'); plot([0 timeWindowInS],[StimPeriodRipsInBin StimPeriodRipsInBin]*length(stimSorted_vector),'b'); 
%plot([0 timeWindowInS],[NoStimPeriodRipsInBin NoStimPeriodRipsInBin]*length(stimSorted_vector),'g');   
title(['Peri-stimulus ripple events (select stimuli: trials: ' num2str(min(trials)) '-' num2str(max(trials)) ' freq: ' num2str(min(freqs)) ' - ' num2str(max(freqs)) ' attenuations: ' num2str(min(attns)) '-' num2str(max(attns))]); axis([0 timeWindowInS 0 20]) 


%% Sixth, Make histogram for peri-stim SWR, as above, but using fake windows in non-stim period, as control.

postStimLag = 45*1500; % how much time to wait after stim ends to start considering "non-stim" timepoints?
preEndLag= 5*1500;
stimNu=floor(900*(StimPeriodRipsInBin/NoStimPeriodRipsInBin)); % how many fake stimuli do you want? (it is 
%normalized so that we loook at the same total # of ripples, since ripple rates are different in the two periods)
stimInds_noStim= randi([max(stimInds)+postStimLag, length(allEEGs(:,10))-preEndLag],[stimNu,1]); %get random indices


% Row indices for freq-attn sorted stimuli. 
periStimRip_sorted_fake=NaN(1+samplesPerWindow/bin,length(stimInds_noStim));

%get the fake histograms - compile into an array.
for m=1:length(stimInds_noStim)
    stimTrigSorted_fake=stimInds_noStim(m);
    stimTrigWindowStartSorted_fake=stimTrigSorted_fake-samplesPerWindow/2;
    stimTrigWindowEndSorted_fake=stimTrigSorted_fake+samplesPerWindow/2;
    stim_edgesSorted_fake=stimTrigWindowStartSorted_fake:bin:stimTrigWindowEndSorted_fake;
    periStimRip_sorted_fake(:,m)=histc(ripInds,stim_edgesSorted_fake);
end

%plotting the new histogram. 
periStimRipSorted_mean_fake=mean(periStimRip_sorted_fake,2);
periStimRipSorted_sum_fake=sum(periStimRip_sorted_fake,2);

figure
bar(0:bin_sec:timeWindowInS,periStimRipSorted_sum_fake,'histc'), hold on, plot(timeWindowInS/2,0,'r'); plot([0 timeWindowInS],[NoStimPeriodRipsInBin NoStimPeriodRipsInBin]*length(stimInds_noStim),'b'); 
%plot([0 timeWindowInS],[NoStimPeriodRipsInBin NoStimPeriodRipsInBin]*length(stimSorted_vector),'g');   
title(['Peri-stimulus ripple events (no-stim period, random timepoints): trials: ' num2str(min(trials)) '-' num2str(max(trials)) ' freq: ' num2str(min(freqs)) ' - ' num2str(max(freqs)) ' attenuations: ' num2str(min(attns)) '-' num2str(max(attns))]); axis([0 timeWindowInS 0 14]) 


%% Seventh, computing peri-stim histogram, but using fake stimulus windows WITHIN the stim period.

stimNu=900; % how many fake stimuli do you want?
stimInds_noStim= randi([min(stimInds) max(stimInds)],[stimNu,1]);

% Row indices for freq-attn sorted stimuli. 
periStimRip_sorted_fake=NaN(1+samplesPerWindow/bin,length(stimInds_noStim));

%Get histograms collected in arrays
for m=1:length(stimInds_noStim)
    stimTrigSorted_fake=stimInds_noStim(m);
    stimTrigWindowStartSorted_fake=stimTrigSorted_fake-samplesPerWindow/2;
    stimTrigWindowEndSorted_fake=stimTrigSorted_fake+samplesPerWindow/2;
    stim_edgesSorted_fake=stimTrigWindowStartSorted_fake:bin:stimTrigWindowEndSorted_fake;
    periStimRip_sorted_fake(:,m)=histc(ripInds,stim_edgesSorted_fake);
end

%plotting the new histogram
periStimRipSorted_mean_fake=mean(periStimRip_sorted_fake,2);
periStimRipSorted_sum_fake=sum(periStimRip_sorted_fake,2);

figure
bar(0:bin_sec:timeWindowInS,periStimRipSorted_sum_fake,'histc'), hold on, plot(timeWindowInS/2,0,'r'); plot([0 timeWindowInS],[StimPeriodRipsInBin StimPeriodRipsInBin]*length(stimInds_noStim),'b'); 
%plot([0 timeWindowInS],[NoStimPeriodRipsInBin NoStimPeriodRipsInBin]*length(stimSorted_vector),'g');   
title(['Peri-stimulus ripple events (no-stim period, random timepoints): trials: ' num2str(min(trials)) '-' num2str(max(trials)) ' freq: ' num2str(min(freqs)) ' - ' num2str(max(freqs)) ' attenuations: ' num2str(min(attns)) '-' num2str(max(attns))]); axis([0 timeWindowInS 0 14]) 



%% PART II: STIM TRIGGERED LFP, AND CORRELATIONS. First, calculate average spectrogram and EEG triggered by all stimuli

%There is commented out script here for calculating spectrograms, but I did
%not complete them.  

% %Spec parameters
% SpecWindow_sec=0.4;
% SpecWindow_samples=1500*SpecWindow_sec;

%LFP parameters
LFPwindow_sec=1;
LFPwindow_samples=LFPwindow_sec*1500;


% Gathering EEGs
for tets=1:21;
    stimTrigLFP= zeros(LFPwindow_samples+1,1);
    stimTrigSpecgram= zeros(201,29);
for n=1:length(stimInds);
stimTrig=stimInds(n);

%SWR centered spectrogram
%WindowStart=stimTrig-SpecWindow_samples/2;
%WindowEnd=stimTrig+SpecWindow_samples/2;
%[S,F,T,P] = spectrogram(allEEGs(WindowStart:WindowEnd,tets), 25,5, 0:200, 1500);
%stimTrigSpecgram=(stimTrigSpecgram+P)/n;

%swr centered lfp
LFPwindowStart=stimTrig-LFPwindow_samples/2;
LFPwindowEnd=stimTrig+LFPwindow_samples/2;
stimTrigLFP=stimTrigLFP+allEEGs([LFPwindowStart:LFPwindowEnd],tets)/n;

end

figure; 
subplot(2,1,2), plot(stimTrigLFP); hold on, plot(zeros(1,3001)+LFPwindow_samples/2,[-1500:1500], 'r--'); axis([0 1600 -max(stimTrigLFP)-100 max(stimTrigLFP)+100]); hold on

end

%% Second, Looking at LFPs only during specific stimuli trials, freqs, and attns that you are interested in.

%Manually input these variables:
tetset_big = [5 7 10 14 16 19] %empirically determined "good tetrodes" for borg day 17 I used 7,10,16 for the three regions.
tetrodes=[1:7 9:11 14 16:20]; % you can use tetset_big if you want.  

%YOU SHOULD LEAVE THESE PARAMETERS AS MIN:MAX, TO GET COMPLETE ARRAY.  HOWEVER, SINCE ALL ATTNS
%WILL BE AVERAGED, SHOULD ONLY CHOOSE THE ATTNS YOU ARE INTERESTED IN.  
trials2=[1:20]
freqs2=[1:15]
attns2=3:5 

% Calculating frequency-atten indices
FrAtInd_list2=[];
clear n; clear m;
for n=1:length(freqs2);
    fr=freqs2(n);
    FrAtInd2=zeros(1,length(attns2));
    for at=attns2
        FrAtInd2=(fr-1)*5+at;
        FrAtInd_list2=[FrAtInd_list2 FrAtInd2];
    end
end

sort(FrAtInd_list2);


LFPwindow_sec=10; %window size in sec
LFPwindow=LFPwindow_sec*1500; %window size in samples - can be big if want to use "control" samples outside of peri-stim area.  see where the array is used below, for correlation measurements.
LFPdurStim= NaN(max(tetrodes),max(trials2), max(FrAtInd_list2),LFPwindow+1);

%Compile all of the windowed EEGs.  
clear FrAt, clear t, clear tet;
for tet=tetrodes;
    for t=trials2;
        for FrAt=FrAtInd_list2;
            windowStim=(allTrialsPerStimInd(t,FrAt)-LFPwindow/2):(allTrialsPerStimInd(t,FrAt)+LFPwindow/2);
            LFPdurStim(tet,t,FrAt,:)=allEEGs(windowStim,tet);
        end
    end
end

LFPwindow_real= 1*1500; %window size to actually consider for these analyses. The one above was bigger just to get all that extra stuff in the arrays, in case you need it.

% new array with each cell corresponding to one tetrode - LFP averaged over all trials, freq, and atten computed above. 
LFPtrial_mean=NaN(max(tetrodes),max(trials2),LFPwindow_real+1);
LFPtet_mean=NaN(max(tetrodes),LFPwindow_real+1);

%averaging EEGs
for tet=tetrodes;
    for t=trials2           
        window2=LFPwindow/2 - LFPwindow_real/2:LFPwindow/2 + LFPwindow_real/2;
           LFPtrial_mean(tet,t,:)= squeeze(mean(LFPdurStim(tet,t,FrAtInd_list2,window2),3));   %averages over all Frequency-attenuation indices - combined out of convenience, since frq and attn are in same dimension
               
    end
    LFPtet_mean(tet,:)=squeeze(mean(LFPtrial_mean(tet,trials2,:),2));% averages over all trials, to get one LFP per tetrode
end

%plot the average EEG responses.  
close all
for tet=tetrodes;
   figure(tet); hold on
   plot(0:(LFPwindow_real/1500)/LFPwindow_real:LFPwindow_real/1500,squeeze(LFPtet_mean(tet,:)),'Color',0.5*[1 1 1]); axis([0 LFPwindow_real/1500 min(LFPtet_mean(10,:))-40 max(LFPtet_mean(10,:))+10]); hold on; plot(LFPwindow_sec/2,0,'r^'), hold on; title({['average peri-stim LFP for tetrode:' num2str(tet) 'over trials:' num2str(min(trials2)) '-' num2str(max(trials2))] ['freqs:' num2str(freqs2)] ['Attns:' num2str(attns2)]}); hold on
end



%% Third, does the EEG response attenuate over time (i.e. over trial #s)?
% looking at trial-by-trial LFP (i.e. no mean as in above) and taking integral of stimulus response.  

LFPdurStim; % Take this array (made above) with EEG trace for each tetrode, trial, frequency, and attenuation

% want to plot trial# against average LFP integral for each frequency (collapse attenuations into one frequency)

LFP_respInt=NaN(max(tetrodes),max(trials2),max(freqs2));
Int_length_sec=0.3; %how much EEG to integrate as your "response" post-stim.  
IntLength=Int_length_sec*1500;

%calculate the integrals
for tet=tetrodes;
    for t=trials2;
        for f=freqs2
            AllAttnPerFreq=f*5-4 +attns2-1; % refers to all attenuation stimuli for a single frequency, f, which allows us to average LFP over all attenuations, for each f.
            LFP_respInt(tet,t,f)=mean(squeeze(mean(abs(LFPdurStim(tet,t,AllAttnPerFreq,LFPwindow/2:(LFPwindow/2+IntLength))),3)));%first averages LFP (use the baseline of "0") overall attenuations for each frequency, then makes everything positive, then takes mean over time(i.e. integral/width)
        end
    end
end
    

%plot for each tetrode: trial vs. average integral of EEG response to stim.
for tet=tetrodes
    figure
    imagesc(squeeze(LFP_respInt(tet,:,:)))
     hold on; colorbar
    title({['Stim response LFP integral for tetrode:' num2str(tet) 'yaxis trials increase going down. trials:' num2str(min(trials2)) '-' num2str(max(trials2))] ['x axis freqs (all atten. averaged):' num2str(freqs2)]})
   
end

%% Fourth: GET AVERAGE RESPONSE INTENSITIES FOR SPECIFIC STIMULI.

tetrodes9=[7 10 16];

freqs9=freqs2 % input here frequencies averaged over. default is all freqs in freqs2
trials9= 11:20; %input here which trials you want

clear LFPIntMean; clear LFPIntSTD;

%take mean and std over all trials.
for tet=tetrodes9
        LFPIntMean(tet)=mean(mean(LFP_respInt(tet,trials9,freqs9)));
        LFPIntSTD(tet)=std(reshape(LFP_respInt(tet,trials9,freqs9),numel(LFP_respInt(tet,trials9,freqs9)), 1));
end

%I can't remember what the below commented out script was for. 
%figure
%bar(tetrodes9,LFPIntMean(LFPIntMean>0)); set(gca,'XTick',3)

LFPIntMean=LFPIntMean(LFPIntMean>0)
LFPIntSTD=LFPIntSTD(LFPIntSTD>0)

%% Fifth, Average the above values over all frequencies for each trial.  Is there a solely trial-dependent effect?
    
%plot response integrals for all stimuli in a 2d heat plot.
    figure
    imagesc(squeeze(mean(LFP_respInt(tetrodes,:,freqs2),3)))
     hold on; colorbar
    title({['Stim response LFP integral for tetrodes(yaxis):' num2str(tetrodes)] ['trials(xaxis), all stim response averaged:' num2str(min(trials2)) '-' num2str(max(trials2))]})

%SAME AS THE ABOVE HEAT PLOT, BUT HERE PLOTTING
% ON SEPARATE PLOTS, SO THAT VALUES SCALE WITHIN THEMSELVES AND DIFFERENCES WITHIN EACH BRAIN REGION, ACROSS TRIALS, ARE THUS EASIER TO SEE.

% AC only
    figure; subplot(3,1,1)
    tet_AC=tetrodes(tetrodes<8);
    imagesc(squeeze(mean(LFP_respInt(tet_AC,:,freqs2),3)),[30 145]) 
     hold on; colorbar
    title({['AC ONLY: Stim response LFP integral for tetrodes(yaxis):' num2str(tetrodes)] ['trials(xaxis), all stim response averaged:' num2str(min(trials2)) '-' num2str(max(trials2))]})

% HC only    
    subplot(3,1,2)
    tet_HC=tetrodes(tetrodes>7 & tetrodes<15);
    imagesc(squeeze(mean(LFP_respInt(tet_HC,:,freqs2),3)),[200 450])
    hold on; colorbar
    title({['HC ONLY: Stim response LFP integral for tetrodes(yaxis):' num2str(tetrodes)] ['trials(xaxis), all stim response averaged:' num2str(min(trials2)) '-' num2str(max(trials2))]})

% PFC only
    subplot(3,1,3)
    tet_PFC=tetrodes(tetrodes>14);
    imagesc(squeeze(mean(LFP_respInt(tet_PFC,:,freqs2),3)),[35 250])
    hold on; colorbar
    title({['PFC ONLY: Stim response LFP integral for tetrodes(yaxis):' num2str(tetrodes)] ['trials(xaxis), all stim response averaged:' num2str(min(trials2)) '-' num2str(max(trials2))]})
   
  
    
%% Sixth, looking at trial by trial correlation between HC, AC, and PFC. 

% 1) for each window (peri-stimulus) calculate 'corr' of EEG between AC-PFC, HC-PFC, and HC-AC (choose one tetrode each).
% 2) Have an imagesc plot for each of those three relationships, 2 dimensions = trial and frequency
% 3) Look at both trial by trial and average over many freqs/attn.

% Input variables
CWindow_sec=0.5; %Stim-triggered window
CorrWindow=CWindow_sec*1500;
offSetSwitch= 0;% Enter 1 to use an offset in the window.  Enter 0 for no offset. see directly below
offSetSize= 5000;%Value to shift window to the left by.  Enter 0 normally, to look at post-stim.  Enter "-CorrWindow" to get pre-stim.  Enter large number (e.g. 5000) to get "control" sample.  Need to make sure that "LFPdurStim" was made using window sizes of required magnitude.
offSetWindow= offSetSwitch*offSetSize; 
AC_tet= 7; % use single tetrodes with strongest LFPs for first pass analysis.
HC_tet= 10;
PFC_tet= 16; 
LFPdurStim; %Rewriting here.  EEGS: (tetrodes, trials, Freq-att index, EEG data) - see above

trials3=1:20;
freqs3=1:15;
attns3=3:5;
AC_HC_PFC_corrCell=cell(max(freqs3), max(attns3));

% To make a cell with corr coeff matrix nxn where n=number of EEG traces. cells are freq, attn.  trials are columns in this matrix, with min(t)-max(t) repeating 3 times total, AC-HC-PFC.  
    for f=freqs3;
        for a=attns3;
            FrAt=f*5-4+a-1; %an index that combines freq and attn info, since they are in same row of array with stimuli indices.
            AC_HC_PFC_corrCell{f,a}= corr([(squeeze(LFPdurStim(AC_tet,trials3,FrAt,LFPwindow/2-offSetWindow:LFPwindow/2+CorrWindow-offSetWindow)))' (squeeze(LFPdurStim(HC_tet,trials3,FrAt,LFPwindow/2-offSetWindow:LFPwindow/2+CorrWindow-offSetWindow)))' (squeeze(LFPdurStim(PFC_tet,trials3,FrAt,LFPwindow/2-offSetWindow:LFPwindow/2+CorrWindow-offSetWindow)))']);
        end
    end
    
    
%Plot, for each trial, correlations between the three pairs of regions.

% AC_HC_PFC_corrCell{f,a} is 2d, comparing correlation for different trials and different brain areas.
% We are only interested in the same trial, across
% brain regions, so here we will isolate just the diagonals of the matrix.  The diagonals
% correspond to the area1-area2 pairs indicated in the names.  The indices
% below pick out the diagonals
AcHc_sameTrialsx=1:20;
AcHc_sameTrialsy=21:40;
AcPfc_sameTrialsx=1:20;
AcPfc_sameTrialsy=41:60;    
HcPfc_sameTrialsx=21:40;
HcPfc_sameTrialsy=41:60;

% stims that you want to plot.  each freq-attn combo has own figure.
freqs4=1:15;
attns4=5;

% Does the presence or absense of a ripple in the stimulus window affect
% the correlation values we get?  Here, make an analogous array to stick in index numbers 
% for ripples, so I can ask whether the presence of ripples correlates with the 3 correlations between 3 pairs of areas
Stim4rips=NaN(max(freqs4),max(attns4),max(trials3));
ripsAtEachStim=NaN(max(freqs4),max(attns4),max(trials3));


%Finally, plot correlation between each of three pairs as bar graphs, and
%plot as a line on the same chart whether a ripple was present.  
for f=freqs4
    for a=attns4;
    Stim4rips(f,a,trials3)=allTrialsPerStimInd(trials3,f*5-5+a);
        for t=trials3
            ripsAtEachStim(f,a,t)=length(find(ripInds<Stim4rips(f,a,t)-offSetWindow & ripInds>(Stim4rips(f,a,t)-2*CorrWindow-offSetWindow))); %getting the numbers of ripples after each trial for specific f and a
        end
        HcPfcCorr(f,a,trials3)=squeeze(diag(AC_HC_PFC_corrCell{f,a}(HcPfc_sameTrialsx,HcPfc_sameTrialsy),0)); % getting area-area correlation for each trial and stim.
        AcHcCorr(f,a,trials3)=squeeze(diag(AC_HC_PFC_corrCell{f,a}(AcHc_sameTrialsx,AcHc_sameTrialsy),0));
        AcPfcCorr(f,a,trials3)=squeeze(diag(AC_HC_PFC_corrCell{f,a}(AcPfc_sameTrialsx,AcPfc_sameTrialsy),0));
        
        figure; bar([trials3' trials3' trials3'],[squeeze(HcPfcCorr(f,a,trials3)) squeeze(AcHcCorr(f,a,trials3)) squeeze(AcPfcCorr(f,a,trials3))]); hold on; plot(trials3, (squeeze(ripsAtEachStim(f,a,trials3))/2),'r'); xlabel('trial #'); ylabel('Correlation Coefficient'); title({['freq:' num2str(f) 'attn:' num2str(a)],'blue=HC-PFC; green=AC-HC; red=AC-PFC'}) ;  hold on; axis tight
    end
end



%% Seventh, find average correlations for three pairs for stim windows with and without ripples (ripples either before or after stim)

%arrays to collect correlation values in separate arrays depending on whether
% a ripple was present during stimulus window.  
                HcPfc_WithRip=NaN(max(freqs4),max(attns4), max(trials3));
                AcHc_WithRip=NaN(max(freqs4),max(attns4), max(trials3));
                AcPfc_WithRip=NaN(max(freqs4),max(attns4), max(trials3));
                HcPfc_WithoutRip=NaN(max(freqs4),max(attns4), max(trials3));
                AcHc_WithoutRip=NaN(max(freqs4),max(attns4), max(trials3));
                AcPfc_WithoutRip=NaN(max(freqs4),max(attns4), max(trials3));

%compile correlation values.
for f=freqs4;
    for a= attns4;
        for t=trials3;
            if ripsAtEachStim(f,a,t)>0;                             
                HcPfc_WithRip(f,a,t)=HcPfcCorr(f,a,t);
                AcHc_WithRip(f,a,t)=AcHcCorr(f,a,t);
                AcPfc_WithRip(f,a,t)=AcPfcCorr(f,a,t);
            else
                HcPfc_WithoutRip(f,a,t)=HcPfcCorr(f,a,t);
                AcHc_WithoutRip(f,a,t)=AcHcCorr(f,a,t);
                AcPfc_WithoutRip(f,a,t)=AcPfcCorr(f,a,t);
            end
        end
    end
end


%find mean correlations for above 6 pairs (i.e. 3, each w. or w/o ripples)

HcPfcMeanCorr_WithRip=mean(HcPfc_WithRip(HcPfc_WithRip>0))
AcHcMeanCorr_WithRip=mean(AcHc_WithRip(AcHc_WithRip>0))
AcPfcMeanCorr_WithRip=mean(AcPfc_WithRip(AcPfc_WithRip>0))  
HcPfcMeanCorr_WithoutRip=mean(HcPfc_WithoutRip(HcPfc_WithoutRip>0))
AcHcMeanCorr_WithoutRip=mean(AcHc_WithoutRip(AcHc_WithoutRip>0))
AcPfcMeanCorr_WithoutRip=mean(AcPfc_WithoutRip(AcPfc_WithoutRip>0))

%SHOULD CALCULATE S.E.M., NOT JUST STD AS I DID HERE.
%SHOULD CALCULATE S.E.M., NOT JUST STD AS I DID HERE.
%SHOULD CALCULATE S.E.M., NOT JUST STD AS I DID HERE.

HcPfcSTDCorr_WithRip=std(HcPfc_WithRip(HcPfc_WithRip>0))
AcHcSTDCorr_WithRip=std(AcHc_WithRip(AcHc_WithRip>0))
AcPfcSTDCorr_WithRip=std(AcPfc_WithRip(AcPfc_WithRip>0))   
HcPfcSTDCorr_WithoutRip=std(HcPfc_WithoutRip(HcPfc_WithoutRip>0))
AcHcSTDCorr_WithoutRip=std(AcHc_WithoutRip(AcHc_WithoutRip>0))
AcPfcSTDCorr_WithoutRip=std(AcPfc_WithoutRip(AcPfc_WithoutRip>0))

%plot those blaues, along with STD
bar([1 2 3],[HcPfcMeanCorr_WithRip HcPfcMeanCorr_WithoutRip; AcHcMeanCorr_WithRip AcHcMeanCorr_WithoutRip; AcPfcMeanCorr_WithRip AcPfcMeanCorr_WithoutRip]); title({'Average correlation between areas, sorted by presence of ripple during post-stim. blue=with rip; red=without rip' '1=HC-PFC, 2=AC-HC, 3=AC-PFC'}); hold on; 


%% Eighth, finding correlations between the 3 area correlations (e.g. AC-PFC).  e.g. when AC-PFC is more correlated, does PFC-HC become less correlated?


HcPfcCorr; % Listing the arrays that tell us the area-area correlation for each trial and stim. (freq,atten,trials)
AcHcCorr;
AcPfcCorr;

% stimuli you wish to average over
freqs5=1:15;
attns5=5;
trials5=1:20;

%(this below, just ignore)
% % THIS BELOW ACTUALLY CALCULATES CORRELATION FOR ONE VECTOR CONTAINING ALL
% % PRESENTATIONS OF THE STIMULUS.  A MIX BETWEEN SIGNAL AND NOISE
% % CORRELATION.  SEE BELOW THIS FOR NOISE CORRELATION 
% Putting all correlation values into one column, one cell for each freq,att,trial. Important that stim order is same for all three vectors.  
% HcPfcCorrLinear = reshape(HcPfcCorr(freqs5,attns5,trials5),numel(HcPfcCorr(freqs5,attns5,trials5)),1);
% AcHcCorrLinear= reshape(AcHcCorr(freqs5,attns5,trials5),numel(AcHcCorr(freqs5,attns5,trials5)),1);
% AcPfcCorrLinear= reshape(AcPfcCorr(freqs5,attns5,trials5),numel(AcPfcCorr(freqs5,attns5,trials5)),1);
% 
% This gets you one 3*3 matrix telling you the correlation coeff for each stimulus presentation (index values) btw 2 brain regions.
% [CorrCorr_AcHcPfc CorrCorr_AcHcPfc_Pvalues] = corr([AcHcCorrLinear AcPfcCorrLinear HcPfcCorrLinear])


%(again, ignore this and use the code directly above)
%ACTUAL NOISE CORRELATION, SEE COMMENT IN CAPS ABOVE. (HERE WE SUBTRACT MEAN AND
%DIVIDE BY STD BEFORE COMPILE EACH FREQ/ATTN(20 TRIAL VALUES) INTO ONE
%GIANT VECTOR FOR ALL STIM PRESENTATIONS).i.e. z-score each stim's 20
%trials before compiling into one vector.

% %zscoring
% HcPfcCorr_zscored=NaN(max(freqs5),max(attns5),max(trials5));
% AcHcCorr_zscored=NaN(max(freqs5),max(attns5),max(trials5));
% AcPfcCorr_zscored=NaN(max(freqs5),max(attns5),max(trials5));
% 
% for f=freqs5;
%     for a= attns5;
%  
%         HcPfcCorr_zscored(f,a,trials5)=zscore(HcPfcCorr(f,a,trials5));
%         AcHcCorr_zscored(f,a,trials5)=zscore(AcHcCorr(f,a,trials5));
%         AcPfcCorr_zscored(f,a,trials5)=zscore(AcPfcCorr(f,a,trials5));
%     end
% end
% 
% HcPfcCorrLinear_zscored = reshape(HcPfcCorr_zscored(freqs5,attns5,trials5),numel(HcPfcCorr_zscored(freqs5,attns5,trials5)),1);
% AcHcCorrLinear_zscored= reshape(AcHcCorr_zscored(freqs5,attns5,trials5),numel(AcHcCorr_zscored(freqs5,attns5,trials5)),1);
% AcPfcCorrLinear_zscored= reshape(AcPfcCorr_zscored(freqs5,attns5,trials5),numel(AcPfcCorr_zscored(freqs5,attns5,trials5)),1);
% 
% % This gets you one 3*3 matrix telling you the correlation coeff for each stimulus presentation (index values) btw 2 brain regions.
% [CorrCorr_AcHcPfc_zscored CorrCorr_AcHcPfc_Pvalues_zscored] = corr([AcHcCorrLinear_zscored AcPfcCorrLinear_zscored HcPfcCorrLinear_zscored])





%DOING NOISE CORRELATION A THIRD WAY: THIS TIME ONLY SUBTRACTING MEAN FROM
%EACH STIM'S 20 TRIALS, AND THEN COMPILING INTO 1D VECTOR (I.E. NOT
%DIVIDING BY STD, LIKE LAST TIME THROUGH ZSCORING.

%subtracting mean from each stim 20 trials.
HcPfcCorr_noise=NaN(max(freqs5),max(attns5),max(trials5));
AcHcCorr_noise=NaN(max(freqs5),max(attns5),max(trials5));
AcPfcCorr_noise=NaN(max(freqs5),max(attns5),max(trials5));

for f=freqs5;
    for a= attns5;
 
        HcPfcCorr_noise(f,a,trials5)=bsxfun(@minus,HcPfcCorr(f,a,trials5),mean(HcPfcCorr(f,a,trials5)));
        AcHcCorr_noise(f,a,trials5)=bsxfun(@minus,AcHcCorr(f,a,trials5),mean(AcHcCorr(f,a,trials5)));
        AcPfcCorr_noise(f,a,trials5)=bsxfun(@minus,AcPfcCorr(f,a,trials5),mean(AcPfcCorr(f,a,trials5)));
    end
end

HcPfcCorrLinear_noise = reshape(HcPfcCorr_noise(freqs5,attns5,trials5),numel(HcPfcCorr_noise(freqs5,attns5,trials5)),1);
AcHcCorrLinear_noise= reshape(AcHcCorr_noise(freqs5,attns5,trials5),numel(AcHcCorr_noise(freqs5,attns5,trials5)),1);
AcPfcCorrLinear_noise= reshape(AcPfcCorr_noise(freqs5,attns5,trials5),numel(AcPfcCorr_noise(freqs5,attns5,trials5)),1);

% This gets you one 3*3 matrix telling you the correlation coeff for each stimulus presentation (index values) btw 2 brain regions.
[CorrCorr_AcHcPfc_noise CorrCorr_AcHcPfc_Pvalues_noise] = corr([AcHcCorrLinear_noise AcPfcCorrLinear_noise HcPfcCorrLinear_noise])





%% Ninth, Calculate correlation between brain regions, this time using the integral of the response, not the EEG

%first take absolute values of all traces.  Then take mean of trace to get
%integral.  Then, for each stim (i.e. all 20 trials) subtract the mean, as
%we are looking for noise correlation. (or z-score, to reduce the greater
%weight that would be put on responses from best frequencies (i.e. greater
%values means greater variance around mean).  then put all integral values
%in a #stim_events x (#tetrodes) vector over all stim events and correlate
%(using "corr") over the three tetrodes.

IntLength2 = 0.3*1500; %window size in samples.
tetrodes2=[7 10 16];
trials6=1:20;
freqs6=1:15;
attns6=attns2; %will not iterate over attns, since all attns have been collapsed already into a single frequency value.  (i.e. the EEGs in LFPdurStim already come averaged over all attns in the "attns2" variable.  
%FrAt is f and a combined into one variable.  not ideal.  

offSetSwitch= 0;% Enter 1 to use an offset in the window.  Enter 0 for no offset. see below
offSetSize= 5000;%Value to shift window to the left by.  Enter 0 normally, to look at post-stim.  Enter "-CorrWindow" to get pre-stim.  Enter large number (e.g. 5000) to get "control" sample.  Need to make sure that "LFPdurStim" was made using window sizes of required magnitude.
offSetWindow= offSetSwitch*offSetSize; 
 


% performing the correlation by first subtracting mean from each stim's 20 trials
% so we are actually correlation fluctuation from the mean. 
LFPdurStim_absInt_noise=NaN(max(tetrodes2),max(trials6),max(freqs6));


for tet=tetrodes2;
    for f=freqs6;
        for a=attns6
        LFPdurStim_absInt_noise(tet,trials6,f,a)= bsxfun(@minus, LFPdurStim_absInt(tet,trials6,f,a), mean(LFPdurStim_absInt(tet,trials6,f,a)));
        end
    end
end
       

%The resulting arrays
LFPdurStim_abs; %EEGs, absolute value, averaged over all attns for each freq.
LFPdurStim_absInt; %Same as above, but all points in each EEG averaged to get integral, for each stimulus.  
LFPdurStim_absInt_noise; %mean subtracted from each stim's 20 trials (mean over those trials)
StimTrigIntegral=[];

for tet=tetrodes2;
    StimTrigIntegral(:,tet)= reshape(LFPdurStim_absInt_noise(tet,trials6,freqs6,attns6),numel(LFPdurStim_absInt_noise(tet,trials6,freqs6,attns6)),1); %putting all integral values into columns in a (#stim) x (max tetrodes) array
end

 
[StimRespCorr StimRespCorrPvalues] =corr(StimTrigIntegral(:,tetrodes2))

clear StimTrigIntegral
 
%REPEAT ABOVE, BUT NO SUBTRACTION OF MEAN FIRST.  
% LFPdurStim_abs=NaN(max(tetrodes2),max(trials6),max(freqs6),max(attns6),IntLength2+1);
% LFPdurStim_absInt=NaN(max(tetrodes2),max(trials6),max(freqs6),max(attns6));
% 
% for tet=tetrodes2;
%     for t = trials6;
%         for f=freqs6
%             for a = attns6
%                    FrAt=f*5-4 +a-1;
%                    LFPdurStim_abs(tet,t,f,a,:)= squeeze(abs(LFPdurStim(tet,t,FrAt,LFPwindow/2-offSetWindow:(LFPwindow/2+IntLength2-offSetWindow)))); %eegs, absolute value, averaged over all attns for each freq.
%                    LFPdurStim_absInt(tet,t,f,a)=mean(LFPdurStim_abs(tet,t,f,a,:));
%             end
%         end
%     end
% end
% 
% 
% %SIGNAL CORRELATION
% for tet=tetrodes2;
%     StimTrigIntegral(:,tet)= reshape(LFPdurStim_absInt(tet,trials6,freqs6,attns6),numel(LFPdurStim_absInt(tet,trials6,freqs6,attns6)),1); %putting all integral values into columns in a (#stim) x (max tetrodes) array
% end
% 
%  
% [StimRespCorr StimRespCorrPvalues] =corr(StimTrigIntegral(:,tetrodes2))


% %REPEAT ABOVE, BUT TRY Z-SCORING FOR EACH STIM, AND NOT JUST SUBTRACTING
% %MEAN - to reduce the greater
% %weight that would be put on responses from best frequencies (i.e. greater
% %values means greater variance around mean).  
% 
% 
% for tet=tetrodes2;
%     for f=freqs6;
%         for a = attns6;
%         LFPdurStim_absInt_zscored(tet,trials6,f,a)= zscore(LFPdurStim_absInt(tet,trials6,f,a));
%         end
%     end
%    
% end
% StimTrigIntegral=[];
% %calculating corr (NOTE: OVERRIDES VALUES FROM ABOVE)
% for tet=tetrodes2;
%     StimTrigIntegral(:,tet)= reshape(LFPdurStim_absInt_zscored(tet,trials6,freqs6,attns6),numel(LFPdurStim_absInt_zscored(tet,trials6,freqs6,attns6)),1); %putting all integral values into columns in a (#stim) x (max tetrodes) array
% end
% 
%   
% [StimRespCorr StimRespCorrPvalues] = corr(StimTrigIntegral(:,tetrodes2))
% clear StimTrigIntegral


%% (IGNORE - LEFTOVER CODE THAT I'M NOT USING.)
% %% Compare the average # of SWR after stim with SWR after random times
%  
% 
% % mean SWR in non-stim equivalent window (shift stimTrig above by random #)
% noStimRip_all = NaN(length(stimInds),2);
% for n=1:1477; %stop at 1477 so that adding random number below never leaves epoch (n=1477 is >30000 less than max index)
% noStimTrig=stimInds(n)+rand(1)*30000; %random shift each trial
% noStimTrigWindowEnd=noStimTrig+samplesPerWindow;
% howManyRipplesAfterNoStim=length(find(ripInds>noStimTrig & ripInds<noStimTrigWindowEnd)); % finds the number of ripple events that had indices between low and high defined for trigger window, puts in array.
% noStimRip_all(n,:)=[noStimTrig howManyRipplesAfterNoStim]; % Gives array with trigger Indices and corresponding # of SWR in post-stim window
% end
% 
% % mean value of SWR, either triggered by stim or random index, over all ind
% 
% noStimRip_mean = mean(noStimRip_all(:,2))
% 
% % plot(noStimRip_all(:,2)); figure
% % plot(stimTrigRip_all(:,2))
% 
% % RESULT: seems to not be any difference even when trying window values of 0.2 to 2, or changing random