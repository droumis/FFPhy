%%  Coherence calculated in non-overlapping window segmentation, state 1 vs. state 2 spectra

%%  ingredients:
%    1. eegstruct
%    2. pos 

animalname = 'Bashir';                              % to label graphs later
windowsize=1;                                     % set to taste
pos_Fs=29.97003;                                    % ETSC video standard
    windowsize_possamp=round(windowsize*pos_Fs);
days = 4;                                           % days to analyze
tetrodes = 1:14;                                    % tetrodes to analyze
groundref_tet=1;
epochs = 2;                                         % epochs to analyze
velocity_state_1 = 2;                               % state 1 (immobile)
velocity_state_2 = 8;                               % state 2 (run)

Fs=eegstruct{days(end)}{epochs(end)}{tetrodes(end)}.samprate;   % probably just 1500
    windowsize_samp=round(windowsize*Fs);

%% First, classify windows into state 1 vs. 2 (and mixed states, 3)

windowstates=cell(days(end),epochs(end));  % epochs
statecounts=zeros(3,epochs(end),days(end));     % 3 states over 7 epochs -- this keeps track of number of windows for each

for d=days      
    for e=epochs                                             % 7 epochs / day
        epochvels=pos{d}{e}.data(:,5);
        nowindows=floor(length(epochvels)/windowsize_possamp);
        windowstates{d,e}=nan(nowindows,1);
        for w=1:nowindows
            meanvel=mean(epochvels((1+(w-1)*windowsize_possamp):(w*windowsize_possamp)));
            if meanvel<velocity_state_1                                               % state 1
                windowstates{d,e}(w)=1;
                statecounts(1,e,d)=statecounts(1,e)+1;
            elseif meanvel>velocity_state_2                                           % state 2
                windowstates{d,e}(w)=2;    
                statecounts(2,e,d)=statecounts(2,e)+1;
            else
                windowstates{d,e}(w)=0;                                                % mixed state
                statecounts(3,e,d)=statecounts(3,e)+1;
            end
        end
    end
end

%% Second, collect raw eegs into windows classified by state.
   %% Critically, DE-REFERENCE all windows by adding ground reference tetrode EEG.

eegs=cell(days(end),epochs(end),tetrodes(end),3);

for d=days
    for e=epochs
        eeg_startindex=lookup(pos{d}{e}.data(1,1),(eegstruct{d}{e}{14}.starttime):1/Fs: ...
                                        (eegstruct{d}{e}{14}.starttime+length(eegstruct{d}{e}{14}.data)/Fs));
        for t=tetrodes
            for w=1:(length(windowstates{d,e})-5)                   %% iterate over each nonoverlapping window
                true_time=(w-1)*windowsize_possamp/pos_Fs;    %% actual time at beginning of window
                start_index=eeg_startindex+round(true_time*Fs);

                groundeeg=eegstruct{d}{e}{groundref_tet}.data(start_index:  ...  
                                                             (start_index+windowsize_samp))';
                tetrodeeeg=eegstruct{d}{e}{t}.data(start_index:  ...  
                                                  (start_index+windowsize_samp))';

                if windowstates{d,e}(w)==1
                    eegs{d,e,t,1}=[eegs{d,e,t,1} ; tetrodeeeg+groundeeg];
                elseif windowstates{d,e}(w)==2
                    eegs{d,e,t,2}=[eegs{d,e,t,2} ; tetrodeeeg+groundeeg];
                else
                    eegs{d,e,t,3}=[eegs{d,e,t,3} ; tetrodeeeg+groundeeg];
                end
                
            end
        end
    end
end


%% Third, pool windows by epoch type (run vs. sleep) or by velocity state.

eegs_statepool=cell(days(end),tetrodes(end),3);             %% discards epoch, pools all within state
eegs_epochpool=cell(days(end),tetrodes(end),2);             %% discards states, pools all within epoch 
eegs_daypool=cell(days(end),tetrodes(end)); 

for d=days
    for t=tetrodes
        
        for s=1:3
            dummy=[];
            for e=epochs
                dummy=cat(1,dummy,eegs{d,e,t,s});
            end
            eegs_statepool{d,t,s}=dummy;
            eegs_daypool{d,t}=cat(1,dummy,eegs_daypool{d,t});
        end
        
        for e=epochs
            dummy=[];
            for s=1:3
                dummy=cat(1,dummy,eegs{d,e,t,s});
            end
            eegs_epochpool{d,e,t}=dummy;
        end
        
    end
end

%% Variable list.

windowstates;
statecounts;
eegs;
eegs_statepool;             
eegs_epochpool;             
eegs_daypool; 


%% calculate coherence

params.trialave=1;          %% average over each 2 second "trial"/segment
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

d=4;

tet8_state1=eegs_statepool{d,8,1}';  %% CA1, near CA2
tet8_state2=eegs_statepool{d,8,2}';

tet9_state1=eegs_statepool{d,9,1}';  %% CA2
tet9_state2=eegs_statepool{d,9,2}';

tet13_state1=eegs_statepool{d,13,1}';   %% CA1
tet13_state2=eegs_statepool{d,13,2}';

tet11_state1=eegs_statepool{d,11,1}';   %% CA1
tet11_state2=eegs_statepool{d,11,2}';

tet14_state1=eegs_statepool{d,14,1}';   %% CA1
tet14_state2=eegs_statepool{d,14,2}';

tet2_state1=eegs_statepool{d,2,1}';   %% ctx
tet2_state2=eegs_statepool{d,2,2}';

tet12_state1=eegs_statepool{d,12,1}';   %% CA3
tet12_state2=eegs_statepool{d,12,2}';

tet10_state1=eegs_statepool{d,10,1}';   %% DG
tet10_state2=eegs_statepool{d,10,2}';

tet4_state1=eegs_statepool{d,4,1}';   %% CC
tet4_state2=eegs_statepool{d,4,2}';

tet3_state1=eegs_statepool{d,3,1}';   %% deep ctx
tet3_state2=eegs_statepool{d,3,2}';


%% superficial ctx coherence

figure
hold on
set(gca,'xscale','log');

reftet=2;
refregion='ctx';
state=2;

[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet8_state2,tet2_state2,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet9_state2,tet2_state2,params);
plot(f,C,'g');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet13_state2,tet2_state2,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet14_state2,tet2_state2,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet11_state2,tet2_state2,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet12_state2,tet2_state2,params);       %%%%%%%
plot(f,C,'r');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet10_state2,tet2_state2,params);       %%%%%%%
plot(f,C,'m');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet4_state2,tet2_state2,params);       %%%%%%%
plot(f,C,'c--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet3_state2,tet2_state2,params);       %%%%%%%
plot(f,C,'c');
title(sprintf('%s state %d coherence in %d s windows, %s, Day %d, reftet=%d',refregion,state,windowsize,animalname,d,reftet),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');


%% CA3 coherence

% state 1

figure
hold on
set(gca,'xscale','log');

reftet=12;
refregion='CA3';
state=1;

[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet8_state1,tet12_state1,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet9_state1,tet12_state1,params);
plot(f,C,'g');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet13_state1,tet12_state1,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet14_state1,tet12_state1,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet11_state1,tet12_state1,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet10_state1,tet12_state1,params);       %%%%%%%
plot(f,C,'m');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet4_state1,tet12_state1,params);       %%%%%%%
plot(f,C,'c--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet3_state1,tet12_state1,params);       %%%%%%%
plot(f,C,'c');
title(sprintf('%s state %d coherence in %d s windows, %s, Day %d, reftet=%d',refregion,state,windowsize,animalname,d,reftet),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');



% state 2

figure
hold on
set(gca,'xscale','log');

reftet=12;
refregion='CA3';
state=2;

[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet8_state2,tet12_state2,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet9_state2,tet12_state2,params);
plot(f,C,'g');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet13_state2,tet12_state2,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet14_state2,tet12_state2,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet11_state2,tet12_state2,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet10_state2,tet12_state2,params);       %%%%%%%
plot(f,C,'m');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet4_state2,tet12_state2,params);       %%%%%%%
plot(f,C,'c--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet3_state2,tet12_state2,params);       %%%%%%%
plot(f,C,'c');
title(sprintf('%s state %d coherence in %d s windows, %s, Day %d, reftet=%d',refregion,state,windowsize,animalname,d,reftet),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');


%% CA1 coherence

% state 1

figure
hold on
set(gca,'xscale','log');

reftet=14;
refregion='CA1';
state=1;

[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet8_state1,tet14_state1,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet9_state1,tet14_state1,params);
plot(f,C,'g');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet13_state1,tet14_state1,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet14_state1,tet14_state1,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet11_state1,tet14_state1,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet10_state1,tet14_state1,params);       %%%%%%%
plot(f,C,'m');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet4_state1,tet14_state1,params);       %%%%%%%
plot(f,C,'c--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet3_state1,tet14_state1,params);       %%%%%%%
plot(f,C,'c');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet12_state1,tet14_state1,params);       %%%%%%%
plot(f,C,'r');
title(sprintf('%s state %d coherence in %d s windows, %s, Day %d, reftet=%d',refregion,state,windowsize,animalname,d,reftet),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');



% state 2

figure
hold on
set(gca,'xscale','log');

reftet=14;
refregion='CA1';
state=2;

[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet8_state2,tet14_state2,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet9_state2,tet14_state2,params);
plot(f,C,'g');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet13_state2,tet14_state2,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet14_state2,tet14_state2,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet11_state2,tet14_state2,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet10_state2,tet14_state2,params);       %%%%%%%
plot(f,C,'m');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet4_state2,tet14_state2,params);       %%%%%%%
plot(f,C,'c--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet3_state2,tet14_state2,params);       %%%%%%%
plot(f,C,'c');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet12_state2,tet14_state2,params);       %%%%%%%
plot(f,C,'r');
title(sprintf('%s state %d coherence in %d s windows, %s, Day %d, reftet=%d',refregion,state,windowsize,animalname,d,reftet),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');


%% CA2 coherence

% state 1

figure
hold on
set(gca,'xscale','log');

reftet=9;
refregion='CA2';
state=1;

[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet8_state1,tet9_state1,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet13_state1,tet9_state1,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet14_state1,tet9_state1,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet11_state1,tet9_state1,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet10_state1,tet9_state1,params);       %%%%%%%
plot(f,C,'m');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet4_state1,tet9_state1,params);       %%%%%%%
plot(f,C,'c--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet3_state1,tet9_state1,params);       %%%%%%%
plot(f,C,'c');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet2_state1,tet9_state1,params);       %%%%%%%
plot(f,C,'c');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet12_state1,tet9_state1,params);       %%%%%%%
plot(f,C,'r');
title(sprintf('%s state %d coherence in %d s windows, %s, Day %d, reftet=%d',refregion,state,windowsize,animalname,d,reftet),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');


% state 2

figure
hold on
set(gca,'xscale','log');

reftet=9;
refregion='CA2';
state=2;

[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet8_state2,tet9_state2,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet13_state2,tet9_state2,params);
plot(f,C,'g--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet14_state2,tet9_state2,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet11_state2,tet9_state2,params);       %%%%%%%
plot(f,C,'k');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet10_state2,tet14_state2,params);       %%%%%%%
plot(f,C,'m');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet4_state2,tet9_state2,params);       %%%%%%%
plot(f,C,'c--');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet3_state2,tet9_state2,params);       %%%%%%%
plot(f,C,'c');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet2_state2,tet9_state2,params);       %%%%%%%
plot(f,C,'c');
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(tet12_state2,tet9_state2,params);       %%%%%%%
plot(f,C,'r');
title(sprintf('%s state %d coherence in %d s windows, %s, Day %d, reftet=%d',refregion,state,windowsize,animalname,d,reftet),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');

     
     
     
     
     
