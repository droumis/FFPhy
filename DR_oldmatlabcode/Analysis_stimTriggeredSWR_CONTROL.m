% LOOKING AT CORRELATIONS BETWEEN PAIRS OF BRAIN REGIONS, BUT AS A CONTROL
% USING RANDOM TIMEPPOINTS, NOT AUDITORY STIMULI.  CODE COPIED OVER FROM
% "ANALYSIS_STIMTRIGGEREDSWR" AND MODIFIED ACCORDINGLY.

% WILL LOOK AT RANDOM TIMEPOINTS DURING THE POST-STIM PERIOD, WHEN NO
% SOUNDS WERE PRESENTED.

%%  GETTING TIMESTAMPS OF "RANDOM" TIMES IN NON-STIM PERIOD AFTER STIMULI WERE PRESENTED AND FINISHED.


postStimLag = 45*1500; % how much time to wait after stim ends to start considering "non-stim" timepoints?
preEndLag= 5*1500;
stimNu=300; % how many fake stimuli do you want?
stimInds_noStim= randi([max(stimInds)+postStimLag, length(allEEGs(:,10))-preEndLag],[stimNu,1]);

% PICKING OUT THE EEGS FOR THOSE NON-STIM TIMEPOINTS.  


tetset_big = [5 7 10 14 16 19]; %empirically determined "good tetrodes"
tetrodes=tetset_big;

LFPwindow_sec=10; %window size in sec
LFPwindow=LFPwindow_sec*1500; %window size in samples - can be big if want to use "control" samples outside of peri-stim area.  see where the array is used below, for correlation measurements.
LFPdurStim_fake= NaN(max(tetrodes),max(stimNu),LFPwindow+1);


for tet=tetrodes;
    for s=1:stimNu;
            windowStim=(stimInds_noStim(s)-LFPwindow/2):(stimInds_noStim(s)+LFPwindow/2);
            LFPdurStim_fake(tet,s,:)=allEEGs(windowStim,tet);
    end
end


LFPwindow_real= 1*1500; %window size to actually consider for these analyses.

%% MODIFY SCRIPT BELOW IF YOU WANT TO SEE THE LFP TRACES THAT YOU ARE PICKING OUT RANDOMLY.  
% new array with each cell corresponding to one tetrode - LFP averaged over all trials, freq, and atten computed above. 
% LFPtrial_mean=NaN(max(tetrodes),max(trials2),LFPwindow_real+1);
% LFPtet_mean=NaN(max(tetrodes),LFPwindow_real+1);
% 
% for tet=tetrodes;
%    for t=trials2           
%        window2=LFPwindow/2 - LFPwindow_real/2:LFPwindow/2 + LFPwindow_real/2;
%           LFPtrial_mean(tet,t,:)= squeeze(mean(LFPdurStim_fake(tet,t,FrAtInd_list2,window2),3));   %averages over all Frequency-attenuation indices - combined out of convenience, since frq and attn are in same dimension
%               
%    end
%    LFPtet_mean(tet,:)=squeeze(mean(LFPtrial_mean(tet,trials2,:),2));% averages over all trials, to get one LFP per tetrode
% end
% 
% for tet=tetrodes;
%   figure
%   plot(0:(LFPwindow_real/1500)/LFPwindow_real:LFPwindow_real/1500,squeeze(LFPtet_mean(tet,:))); axis tight; hold on; plot(LFPwindow_sec/2,0,'r^'), hold on; title({['average peri-stim LFP for tetrode:' num2str(tet) 'over trials:' num2str(min(trials2)) '-' num2str(max(trials2))] ['freqs:' num2str(freqs2)] ['Attns:' num2str(attns2)]}); hold on
% end



%% Looking at trial by trial correlation between different HC, AC, and PFC. 

% 1) for each window (peri-stimulus) calculate 'corr' between AC-PFC, HC-PFC, and HC-AC (choose one tetrode each).
% 2) Have an imagesc plot for each of those three relationships, 2 dimensions = trial and frequency
% 3) Look at both trial by trial and average over many freqs/attn.

% Input variables


CWindow_sec=0.5; %Stim-triggered window
CorrWindow=CWindow_sec*1500;
offSetSwitch= 0;% Enter 1 to use an offset in the window.  Enter 0 for no offset.
offSetSize= 5000;%Value to shift window to the left by.  Enter 0 normally, to look at post-stim.  Enter "-CorrWindow" to get pre-stim.  Enter large number (e.g. 5000) to get "control" sample.  Need to make sure that "LFPdurStim" was made using window sizes of required magnitude.
offSetWindow= offSetSwitch*offSetSize; 
AC_tet= 7;
HC_tet= 10;
PFC_tet= 16; % use single tetrodes with strongest seeming LFPs for first pass analysis.
LFPdurStim_fake; %EEGS: 

% select values below: NOTE VALUES HAVE TO HAVE BEEN CHOSEN ABOVE ALREADY (I.E. #2 VALUES BELOW).  THE ARRAY WILL NOT CONTAIN DATA OTHER VALUES.
%trials2
%freqs2
%attns2

AC_HC_PFC_corrCell_fake=NaN(3*max(stimNu));

% To make a cell with corr coeff matrix nxn where n=number of EEG traces. cells are freq, attn.  trials are columns in this matrix, with min(t)-max(t) repeating 3 times total, AC-HC-PFC.  


       
AC_HC_PFC_corrCell_fake= corr([(squeeze(LFPdurStim_fake(AC_tet,1:stimNu,LFPwindow/2-offSetWindow:LFPwindow/2+CorrWindow-offSetWindow)))' (squeeze(LFPdurStim_fake(HC_tet,1:stimNu,LFPwindow/2-offSetWindow:LFPwindow/2+CorrWindow-offSetWindow)))' (squeeze(LFPdurStim_fake(PFC_tet,1:stimNu,LFPwindow/2-offSetWindow:LFPwindow/2+CorrWindow-offSetWindow)))']);

    
    

%% Plot, for each trial, correlations between the three pairs of regions.
close all
%indices in AC_HC_PFC_corrCell{f,a} that correspond to the area1-area2 pairs indicated in the names.

AcHc_sameTrialsx=1:stimNu;
AcHc_sameTrialsy=stimNu+1:2*stimNu;
AcPfc_sameTrialsx=1:stimNu;
AcPfc_sameTrialsy=2*stimNu+1:3*stimNu;    
HcPfc_sameTrialsx=stimNu+1:2*stimNu;
HcPfc_sameTrialsy=2*stimNu+1:3*stimNu;


% stims that you want to plot. 

%array to stick in index numbers for ripples, so I can ask whether the presence of ripples correlates with the 3 correlations between 3 pairs of areas
Stim4rips_fake=NaN(max(stimNu),1);
ripsAtEachStim_fake=NaN(max(stimNu),1);



%plot correlation between each of three pairs.  
for s=1:stimNu;
        ripsAtEachStim_fake(s)=length(find(ripInds<stimInds_noStim(s)-offSetWindow & ripInds>(stimInds_noStim(s)-2*CorrWindow-offSetWindow))); %getting the numbers of ripples after each trial for specific f and a
end

        HcPfcCorr_fake=squeeze(diag(AC_HC_PFC_corrCell_fake(HcPfc_sameTrialsx,HcPfc_sameTrialsy),0)); % getting area-area correlation for each trial and stim.
        AcHcCorr_fake=squeeze(diag(AC_HC_PFC_corrCell_fake(AcHc_sameTrialsx,AcHc_sameTrialsy),0));
        AcPfcCorr_fake=squeeze(diag(AC_HC_PFC_corrCell_fake(AcPfc_sameTrialsx,AcPfc_sameTrialsy),0));
   

        %figure; plot(1:stimNu,squeeze(HcPfcCorr_fake),'b',1:stimNu,squeeze(AcHcCorr_fake),'g',1:stimNu, squeeze(AcPfcCorr_fake),'r'); hold on; title({'Random non-stim period, EEG cross-correlation', 'blue=HC-PFC; green=AC-HC; red=AC-PFC'}) ; plot(barstimNu,ripsAtEachStim_fake,'k'), hold on;

                figure; bar([(1:40)' (1:40)' (1:40)'],[squeeze(HcPfcCorr_fake(1:40)) squeeze(AcHcCorr_fake(1:40)) squeeze(AcPfcCorr_fake(1:40))]); hold on; plot(1:40, (squeeze(ripsAtEachStim_fake(1:40))/2),'r'); xlabel('trial #'); ylabel('Correlation Coefficient'); title({['freq:' num2str(f) 'attn:' num2str(a)],'blue=HC-PFC; green=AC-HC; red=AC-PFC'}) ;  hold on; axis tight

        %%
% find average correlation for three pairs for stim windows with and without ripples (either before or after)

%arrays for time stamps with ripples, or without ripples
                HcPfc_WithRip_fake=NaN(stimNu,1);
                AcHc_WithRip_fake=NaN(stimNu,1);
                AcPfc_WithRip_fake=NaN(stimNu,1);
                HcPfc_WithoutRip_fake=NaN(stimNu,1);
                AcHc_WithoutRip_fake=NaN(stimNu,1);
                AcPfc_WithoutRip_fake=NaN(stimNu,1);

for s=1:stimNu;
            if ripsAtEachStim_fake(s)>0;                             
                HcPfc_WithRip_fake(s)=HcPfcCorr_fake(s);
                AcHc_WithRip_fake(s)=AcHcCorr_fake(s);
                AcPfc_WithRip_fake(s)=AcPfcCorr_fake(s);
            else
                HcPfc_WithoutRip_fake(s)=HcPfcCorr_fake(s);
                AcHc_WithoutRip_fake(s)=AcHcCorr_fake(s);
                AcPfc_WithoutRip_fake(s)=AcPfcCorr_fake(s);
            end
        end



% find mean correlations for above 6 pairs (i.e. 3, each w. or w/o ripples)

HcPfcMeanCorr_WithRip_fake=mean(HcPfc_WithRip_fake(HcPfc_WithRip_fake>0))
AcHcMeanCorr_WithRip_fake=mean(AcHc_WithRip_fake(AcHc_WithRip_fake>0))
AcPfcMeanCorr_WithRip_fake=mean(AcPfc_WithRip_fake(AcPfc_WithRip_fake>0))   
HcPfcMeanCorr_WithoutRip_fake=mean(HcPfc_WithoutRip_fake(HcPfc_WithoutRip_fake>0))
AcHcMeanCorr_WithoutRip_fake=mean(AcHc_WithoutRip_fake(AcHc_WithoutRip_fake>0))
AcPfcMeanCorr_WithoutRip_fake=mean(AcPfc_WithoutRip_fake(AcPfc_WithoutRip_fake>0))

%find standard deviations;
HcPfcSTDCorr_WithRip_fake=std(HcPfc_WithRip_fake(HcPfc_WithRip_fake>0))
AcHcSTDCorr_WithRip_fake=std(AcHc_WithRip_fake(AcHc_WithRip_fake>0))
AcPfcSTDCorr_WithRip_fake=std(AcPfc_WithRip_fake(AcPfc_WithRip_fake>0))   
HcPfcSTDCorr_WithoutRip_fake=std(HcPfc_WithoutRip_fake(HcPfc_WithoutRip_fake>0))
AcHcSTDCorr_WithoutRip_fake=std(AcHc_WithoutRip_fake(AcHc_WithoutRip_fake>0))
AcPfcSTDCorr_WithoutRip_fake=std(AcPfc_WithoutRip_fake(AcPfc_WithoutRip_fake>0))

figure; bar([1 2 3],[HcPfcMeanCorr_WithRip_fake HcPfcMeanCorr_WithoutRip_fake; AcHcMeanCorr_WithRip_fake AcHcMeanCorr_WithoutRip_fake; AcPfcMeanCorr_WithRip_fake AcPfcMeanCorr_WithoutRip_fake]); title({'Average correlation between areas, sorted by presence of ripple during post-stim. blue=with rip; red=without rip' '1=HC-PFC, 2=AC-HC, 3=AC-PFC'}); hold on; 
%errorbar([1 2 3], [HcPfcMeanCorr_WithRip_fake HcPfcMeanCorr_WithoutRip_fake; AcHcMeanCorr_WithRip_fake AcHcMeanCorr_WithoutRip_fake; AcPfcMeanCorr_WithRip_fake AcPfcMeanCorr_WithoutRip_fake], [0 1; 0 1; 0 1]);

%% finding correlations, across all stim (freq,attn, trials) between the 3 area correlations (e.g. AC-PFC).  i.e. when AC-PFC is more correlated, does PFC-HC become less?


HcPfcCorr_fake; % Listing the arrays that tell us the area-area correlation for each trial and stim. (freq,atten,trials)
AcHcCorr_fake;
AcPfcCorr_fake;


%SIGNAL CORRELATION

%THIS BELOW ACTUALLY CALCULATES CORRELATION FOR ONE VECTOR CONTAINING ALL
%PRESENTATIONS OF THE STIMULUS.  A MIX BETWEEN SIGNAL AND NOISE
%CORRELATION.  SEE BELOW THIS FOR NOISE CORRELATION 
%Putting all correlation values into one column, one cell for each freq,att,trial. Important that stim order is same for all three vectors.  
HcPfcCorrLinear_fake = reshape(HcPfcCorr_fake,numel(HcPfcCorr_fake),1);
AcHcCorrLinear_fake= reshape(AcHcCorr_fake,numel(AcHcCorr_fake),1);
AcPfcCorrLinear_fake= reshape(AcPfcCorr_fake,numel(AcPfcCorr_fake),1);

% This gets you one 3*3 matrix telling you the correlation coeff for each stimulus presentation (index values) btw 2 brain regions.
[CorrCorr_AcHcPfc__signal_fake CorrCorr_AcHcPfc_Pvalues_signal_fake] = corr([AcHcCorrLinear_fake AcPfcCorrLinear_fake HcPfcCorrLinear_fake])



%ACTUAL NOISE CORRELATION, SEE COMMENT IN CAPS ABOVE. (HERE WE SUBTRACT MEAN AND
%DIVIDE BY STD BEFORE COMPILE EACH FREQ/ATTN(20 TRIAL VALUES) INTO ONE
%GIANT VECTOR FOR ALL STIM PRESENTATIONS).i.e. z-score each stim's 20
%trials before compiling into one vector.

%zscoring
HcPfcCorr_zscored_fake=NaN(200,1);
AcHcCorr_zscored_fake=NaN(200,1);
AcPfcCorr_zscored_fake=NaN(200,1);

 
        HcPfcCorr_zscored_fake=zscore(HcPfcCorr_fake);
        AcHcCorr_zscored_fake=zscore(AcHcCorr_fake);
        AcPfcCorr_zscored_fake=zscore(AcPfcCorr_fake);
   

HcPfcCorrLinear_zscored_fake = reshape(HcPfcCorr_zscored_fake,numel(HcPfcCorr_zscored_fake),1);
AcHcCorrLinear_zscored_fake= reshape(AcHcCorr_zscored_fake,numel(AcHcCorr_zscored_fake),1);
AcPfcCorrLinear_zscored_fake= reshape(AcPfcCorr_zscored_fake,numel(AcPfcCorr_zscored_fake),1);

% This gets you one 3*3 matrix telling you the correlation coeff for each stimulus presentation (index values) btw 2 brain regions.
[CorrCorr_AcHcPfc_zscored_fake CorrCorr_AcHcPfc_Pvalues_zscored_fake] = corr([AcHcCorrLinear_zscored_fake AcPfcCorrLinear_zscored_fake HcPfcCorrLinear_zscored_fake])





%DOING NOISE CORRELATION ANOTHER WAY: THIS TIME ONLY SUBTRACTING MEAN FROM
%EACH STIM'S 20 TRIALS, AND THEN COMPILING INTO 1D VECTOR (I.E. NOT
%DIVIDING BY STD, LIKE LAST TIME THROUGH ZSCORING.

%subtracting mean from each stim 20 trials.
HcPfcCorr_noise_fake=NaN(stimNu,1);
AcHcCorr_noise_fake=NaN(stimNu,1);
AcPfcCorr_noise_fake=NaN(stimNu,1);


        HcPfcCorr_noise_fake=bsxfun(@minus,HcPfcCorr_fake,mean(HcPfcCorr_fake));
        AcHcCorr_noise_fake=bsxfun(@minus,AcHcCorr_fake,mean(AcHcCorr_fake));
        AcPfcCorr_noise_fake=bsxfun(@minus,AcPfcCorr_fake,mean(AcPfcCorr_fake));


HcPfcCorrLinear_noise_fake = reshape(HcPfcCorr_noise_fake,numel(HcPfcCorr_noise_fake),1);
AcHcCorrLinear_noise_fake= reshape(AcHcCorr_noise_fake,numel(AcHcCorr_noise_fake),1);
AcPfcCorrLinear_noise_fake= reshape(AcPfcCorr_noise_fake,numel(AcPfcCorr_noise_fake),1);

% This gets you one 3*3 matrix telling you the correlation coeff for each stimulus presentation (index values) btw 2 brain regions.
[CorrCorr_AcHcPfc_noise_fake CorrCorr_AcHcPfc_Pvalues_noise_fake] = corr([AcHcCorrLinear_noise_fake AcPfcCorrLinear_noise_fake HcPfcCorrLinear_noise_fake])





%% NOISE CORRELATION Calculate correlation between brain regions, for the integral of the response.  

%first take absolute values of all traces.  Then take mean of trace to get
%integral.  Then, for each stim (i.e. all 20 trials) subtract the mean, as
%we are looking for noise correlation. (or z-score, to reduce the greater
%weight that would be put on responses from best frequencies (i.e. greater
%values means greater variance around mean).  then put all integral values
%in a #stim_events x (#tetrodes) vector over all stim events and correlate
%(using "corr") over the three tetrodes.

IntLength2 = 0.3*1500; %window size in samples.
tetrodes2=[7 10 16];

offSetSwitch= 0;% Enter 1 to use an offset in the window.  Enter 0 for no offset.
offSetSize= 5000;%Value to shift window to the left by.  Enter 0 normally, to look at post-stim.  Enter "-CorrWindow" to get pre-stim.  Enter large number (e.g. 5000) to get "control" sample.  Need to make sure that "LFPdurStim" was made using window sizes of required magnitude.
offSetWindow= offSetSwitch*offSetSize; 

LFPdurStim_abs_fake=NaN(max(tetrodes2),max(stimNu),IntLength2+1);
LFPdurStim_absInt_fake=NaN(max(tetrodes2),max(stimNu));

for tet=tetrodes2;
    for s=1:stimNu
                   LFPdurStim_abs_fake(tet,s,:)= squeeze(abs(LFPdurStim_fake(tet,s,LFPwindow/2-offSetWindow:(LFPwindow/2+IntLength2-offSetWindow)))); %eegs, absolute value, averaged over all attns for each freq.
                   LFPdurStim_absInt_fake(tet,s)=mean(LFPdurStim_abs_fake(tet,s,:));
        end
end
%%
% %subtracting mean from each stim's 20 trials
% 
% LFPdurStim_absInt_noise=NaN(max(tetrodes2),max(trials6),max(freqs6));
% 
% for tet=tetrodes2;
%     for f=freqs6;
%         LFPdurStim_absInt_noise(tet,trials6,f)= bsxfun(@minus, LFPdurStim_absInt_fake(tet,trials6,f), mean(LFPdurStim_absInt_fake(tet,trials6,f)));
%     end
% end
%        
% 
% %The resulting arrays
% LFPdurStim_abs_fake; %EEGs, absolute value, averaged over all attns for each freq.
% LFPdurStim_absInt_fake; %Same as above, but all points in each EEG averaged to get integral, for each stimulus.  
% LFPdurStim_absInt_noise; %mean subtracted from each stim's 20 trials (mean over those trials)
% 
% for tet=tetrodes2;
%     StimTrigIntegral_fake(:,tet)= reshape(LFPdurStim_absInt_noise(tet,trials6,freqs6),numel(LFPdurStim_absInt_noise(tet,trials6,freqs6)),1); %putting all integral values into columns in a (#stim) x (max tetrodes) array
% end
% 
%   
% [StimRespCorr_fake StimRespCorrPvalues_fake] =corr(StimTrigIntegral_fake(:,tetrodes2))


%%
%REPEAT ABOVE, BUT TRY Z-SCORING FOR EACH STIM, AND NOT JUST SUBTRACTING
%MEAN - to reduce the greater
%weight that would be put on responses from best frequencies (i.e. greater
%values means greater variance around mean).  

StimTrigIntegral_fake=[];
%calculating corr (NOTE: OVERRIDES VALUES FROM ABOVE)
for tet=tetrodes2;
    StimTrigIntegral_fake(:,tet)= reshape(LFPdurStim_absInt_fake(tet,:),numel(LFPdurStim_absInt_fake(tet,:)),1); %putting all integral values into columns in a (#stim) x (max tetrodes) array
end

 
[StimRespCorr_fake StimRespCorrPvalues_fake] = corr(StimTrigIntegral_fake(:,tetrodes2))





