function [spikes,aspkfile]=sss_spkcut_2ch_basic(files1, files2, sdcuts, stimcode, params, consweep, sweepall, sweepstim)

% sss_autospikecut_1ms_basic({'st1e1_r1';'st1e1_r2}';{'st1e2_r1;'st1e2_r2'},
% params, [0 1 0 1 100], 1, [1:1090], {[1:200 601:800]; [201:600 801:1000]; [1001:1090]});

%%% READ THIS SECTION "CAREFULLY!" FOR FIRST TIME USERS

%%% Modified: Shantanu Jadhav, Nov21, 2005

% SSS_AUTOSPIKECUT_1MS_BASIC is a function to be used as a callback in GSORTM,
% and is written to automatically detect and cut spikes at various thresholds.
% Threshold calculation is in respect to the background noise level before the
% stimulus onset.  The user specifies one standard deviation value and the program
% calculates the mean voltage values and standard deviations.

% THIS FILE CUTS 1.5 MS SPIKES FOR SINGLE CHANNEL ONLY (48 samples at 32kHz)
% MAKING THIS A GENERAL FILE FOR ALL INPUT COMBINATIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MANDATORY %%%%%%%%%%%%%%%%%%%%%%%%%%%

%   FILES is a cell array with the .swp file names written individually in seperate
%   lines. ex. {'r49g5P1'; 'r49g5R2';'r49g5SR1'}.  PLEASE NOTE THAT THE FIRST FILE
%   SHOULD BE THE ONE WHERE THE THRESHOLD INFORMATION WILL BE CALCULATED.

%%%%%%%%%%%%%% OPTIONAL, IF NO INPUT, DEFAULTS WILL BE USED:%%%%%%%%%%%%%%
%%%% FIELDS necessary for SSS_STIM_ARCHIVE are PARAMS, STIMCODE(SPIKES.STIMULUS),
%%%% SWEEPALL(SPIKES.SWEEPALL), SWEEPSTIM(SPIKES.SWEEPSTIM) %%%%%%%%%%%%%%

%   SDCUTS is the values to be used to calculate the thresholds.  For example if as
%   SDCUTS there is only one value, say 3, then the program will calculate the mean
%   plus the three sds as the threshold. Example entries: 3|[3:5]|[3 5 9].
%   USE ONLY ONE VALUE
%
%   STIMCODE: (Tricky) For a super-structure stimulus code, which is higher in
%   hierarchy than igorstim. You need one stimcode for each .swp file input. The
%   default stimcodes are 0,1,2...length(files)-1. If you need multiple
%   stimcodes in one file, NEED TO MANUALLY PUT SWEEPALL AND SWEEPSTIM
%   In future, will have an IGORFIELD for this FOR EACH TRIAL, AND THAT WILL BE
%   DEFAULT. THERE WILL HAVE TO CHANGES TO DEFAULT CALCS OF SWEEPALL AND SWEEPSTIM
%   THEN.
%
%   PARAMS: stimonset, bckwindow, sweepd
%       STIMONSET is the onset of the whisker deflection in ms resolution.
%       BCKWINDOW is the window size where the average voltage and the standard
%       deviation of the voltage will be calculated. Resolution in ms.
%       SWEEPD is the sweep duration in ms resolution.
%       params can be empty

%   CONSWEEP: Default=1 implies sweeps will be pushed for multiple files (almost
%   always the case). You want a unique sweep number for all sweeps (for sorting
%   and isi reasons), each with a igorstim (which is in the SPIKES.ALLIGORSTIM field)
%
%   SWEEPALL: Range of sweeps, 1st to last across all files. Default should
%   be done auto by adding sweeps from files by concatanating them
%
%   SWEEPSTIM: Sweep range for each stimulus: one file => one stimulus =>
%   one sweep range. So auto should be picked easily. Need manual input for
%   1 file: multiple stim, OR IF YOU WANT TO DEFINE THE SWEEP RANGE
%   OR if you want sweep stim code irrespective of spike in sweep or not

%%%% OUTPUT
%%%%   For all spikesin all files, stored in all_firstfilename_thr_sdcuts
%%%%   For each file, indiv file op stored in filenmae_thr_sdcuts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% IMP CHANGES: NOV 20, 2005
%% 1. TRYING TO MAKE THIS A GENERAL FILE FOR ALL COMBINATIONS: A) MULTIPLE
%% FILES WITH STIMCODE FOR EACH FILE (MANUAL OR AUTO), OR B)ONE FILE, 1
%% STIMCODE: EASIEST CASE, ENTER ONLY FILENAME AND SDCUTS, OR
%% C)ONE FILE:MULTIPLE STIMCODES: TRICKEIST CASE: HAVE TO MANUALLY ENTER SWEEPALL
%% AND SWEEPSTIM

%   In future, will have an IGORFIELD for STIMCODE FOR ALL SWEEPS IN A TRIAL, AND
%   THAT WILL BE  DEFAULT. THERE WILL HAVE TO CHANGES TO DEFAULT CALCS OF SWEEPALL AND SWEEPSTIM
%   THEN. JUST USE STIMCODE FOR ALL sweeps in a trial: easy.

% shantanu@ucsd.edu
% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off;

if (nargin < 3), sdcuts=5; end
if (nargin < 4), stimcode=(1:size(files1,1))-1; end %% default stimcode = 0,1,2...lth(files)-1
if (nargin < 5), 
    stimonset=200; bckwindow=100; sweepd=700;
else
    if (isempty(params.stimonset)), stimonset = 200; else, stimonset=params.stimonset; end
    if (isempty(params.bckwindow)), bckwindow = 100; else, bckwindow=params.bckwindow; end
    if (isempty(params.sweepd)), sweepd = 700; else, sweepd=params.sweepd; end
end
if (nargin < 6), consweep=1; end

if (exist('sweepall')), 
    if (isempty(sweepall)),autosweepall=1; ;sweepall=[];else,autosweepall=0;end
else
    autosweepall=1; ;sweepall=[];
end

if (exist('sweepstim')),
    if (isempty(sweepstim)),autosweepstim=1; ;sweepstim=[];else,autosweepstim=0;end
else
    autosweepstim=1; sweepstim=[];
end

%if (exist('sdcuts')==0), sdcuts=5; end
if length(stimcode) < size(files1,1)
    stimcode(end+1:size(files1,1))=0;
elseif length(stimcode) > size(files1,1)
    stimcode(size(files1,1):end)=[];
end

unqstim=unique(stimcode);
% Create variables to use later in the code
longsp=0; bsize=1; spikestarts=14;
sfreq=32; % sampling frequency in kHz.

recname = strtok(files1{1},'.');
% [base1,base2] = strtok(files{1},'L');
% recname=[base1 base2(1:end-1)]; %    Rec Site Name

addsweep=0; addtrial=0;
allspikes.waveforms=[]; allspikes.waveforms_ch1=[];allspikes.waveforms_ch2=[];
allspikes.spiketimes=[]; allspikes.swtimes=[];
allspikes.fstimes=[];allspikes.ftimes=[];
allspikes.sweep = [];allspikes.trial = [];
allspikes.stimulus = []; allspikes.igorstim = [];allspikes.nsweeps = 0;
Alligorstim=[]; Allsweepidx=[];

for lp_f=1:size(files1,1)

    %Ch1
    name1 = strtok(files1{lp_f},'.');
    data=load ([files1{lp_f}]);
    h = waitbar(0,'Reorganizing the data into matrix...');
    data = [-9999 data' -9999];
    cuts = find (data == -9999);
    datmat1 = zeros (length(cuts)-1,(cuts(2)-cuts(1))+1); 

    for lp_c = 1: (length (cuts)-1)
        datmat1 (lp_c,:)= data (cuts(lp_c):cuts(lp_c+1)); 
        waitbar(lp_c/(length (cuts)-1));
        a=1;
    end
    clear data ndata h lp_c
    
    %Ch2
    name2 = strtok(files2{lp_f},'.');
    data=load ([files2{lp_f}]);
    h = waitbar(0,'Reorganizing the data: Ch1 into matrix...');
    data = [-9999 data' -9999];
%     cuts = find (data == -9999);          % Cuts should be same as Channel1
    datmat2 = zeros (length(cuts)-1,(cuts(2)-cuts(1))+1); 
    
    for lp_c = 1: (length (cuts)-1)
        datmat2 (lp_c,:)= data (cuts(lp_c):cuts(lp_c+1)); 
        waitbar(lp_c/(length (cuts)-1));
    end
    clear data ndata h lp_c
    

    % KEY FOR THE DATMAT
    %   (which has as many rows as the neumber of sweeps and as many columns as the (cuts(2)-cuts(1))+1))
    % 1st column    -9999
    % 2nd           Rate
    % 3rd           Number of samples
    % 4th           Sweep no
    % 5th           Trial no
    % 6th           Repetition no
    % 7th           Stimulus number
    % 8:N           Sweep data -raw-
    % N+1           -9999 Seperator

    clear data ndata h lp_c

    % CAN CALCULATE THRESHOLD FOR ALL FILES SEPARATELY: KEEP FOR NOW
    % CALCULATE THRESHOLD CROSSINGS ON 1 CHANNEL ONLY
    %% KEEPING SAME THRESHOLD THROUGHOUT, AS IT IS THE SAME CHANNEL:
    %% DIFFERENT STIM CONDITIONS
    
    if lp_f==1 % if the first cycle then calculate the threshold
        tempdat=datmat1(:,8+(stimonset*sfreq-bckwindow*sfreq):8+(stimonset*sfreq));
        tempthr= [mean(mean(tempdat)) mean(std(tempdat))];

        for lpsdc=1:length(sdcuts)
            thrs(lpsdc)= tempthr(1)+sdcuts(lpsdc)*tempthr(2);
        end
    end

    % thrs = 9;

    for lpthrs=1:length(thrs);
        thr= thrs(lpthrs);
        datmat=[]; fstimes=[]; shape1=[]; indct=1; fwaves1=[];
        datmat22=[]; fstimes2=[]; shape2=[]; indct2=1; fwaves2=[];
        h = waitbar(0,'Creating the waveforms...');

        %keyboard
        % Time to find the waveforms and their corresponding spike times.
        % First create a for loop to go over the sweep data (sw).
        for lp_sw = 1: size (datmat1,1)

            %Create a temporary vector where you can write a single sweep data
            tempvec = datmat1(lp_sw,:);
            temp = tempvec (8:end-1); % clean the non A/D values in the temp
            
            %Ch2  
            tempvec2 = datmat2(lp_sw,:);
            temp2 = tempvec2 (8:end-1);
            
            %% SAVE IGORSTIM NUMBER IRRESPECTIVE OF SPIKE OR NOT %%%
            Alligorstim = [Alligorstim tempvec(1,7)];
            Allsweepidx = [Allsweepidx tempvec(1,4)+1];  %%% Save sweep no irrespective
                                                         %%% of spike or not in sw
                                                         %%% For sweepall & sweepstim
            
            % Go and find the threshold passing events on Ch1

            if length(find(temp>=thr))>=1
                g_sp= find (temp>=thr); % values over threshold
                g2_sp= g_sp([find(diff(g_sp)>1)+1]); % threshold crossings in +/- direction
                g2_sp=[g_sp(1) g2_sp]; %spikes start

                g2_sp=g2_sp(find ((g2_sp>spikestarts+1) & (g2_sp<size(tempvec,2)-(56-spikestarts))));
                % get those spikes that are entirely within the 1.5 ms window
                g2_spt= round(g2_sp * 0.3125);    % 0.3125 constant converts the spiketimes from 32kHz sampling to 0.1 ms

                %stc(lp_sw)=length(g2_spt);
                % 1.5 ms spikes
                for lp_sp=1:length(g2_sp);
                    shape1= temp(g2_sp(lp_sp)-spikestarts:g2_sp(lp_sp)+(47-spikestarts));
                    fwaves1= [fwaves1; shape1];
                    fspiketime=(lp_sw*10000)+g2_spt(lp_sp); fstimes=[fstimes; fspiketime];
                    datmat(indct,1)=-9999; % seperator
                    datmat(indct,2)=tempvec(1,4) + 1; % sweep number
                    datmat(indct,3)=tempvec(1,5); % trial number
                    datmat(indct,4)=tempvec(1,6); % repetition number
                    datmat(indct,5)=tempvec(1,7); % stimulus number/igorstim
                    datmat(indct,6)=g2_spt(lp_sp); % time stamp
                    indct=indct+1;
                end
                
                % For Ch2
                for lp_sp=1:length(g2_sp);
                    shape2= temp2(g2_sp(lp_sp)-spikestarts:g2_sp(lp_sp)+(47-spikestarts));            
                    fwaves2= [fwaves2; shape2];
%                     fspiketime2=(lp_sw*10000)+g2_spt(lp_sp); fstimes2=[fstimes2; fspiketime2];  
%                     datmat22(indct2,1)=-9999; % seperator 
%                     datmat22(indct2,2)=tempvec2(1,4) + 1; % sweep number
%                     datmat22(indct2,3)=tempvec2(1,5); % trial number
%                     datmat22(indct2,4)=tempvec2(1,6); % repetition number
%                     datmat22(indct2,5)=tempvec2(1,7); % stimulus number
%                     datmat22(indct2,6)=g2_spt(lp_sp); % time stamp
%                     indct2=indct2+1;
                end
                
            end %end if
            waitbar(lp_sw/size (datmat2,1));
        end
        
        datmat= [datmat fwaves1 datmat(:,1)]; ftimes= [datmat(:,6)];
%       datmat22= [datmat22 fwaves22 datmat22(:,1)]; ftimes22= [datmat22(:,6)];
        
        %scthreshold=sdcuts(lpthrs);
        %fstudied= files{lp_f};

        nsweeps = length (cuts)-1;
        ndatmat=datmat;
        name = strtok(files1{1},'_');

        spikes.waveforms_ch1=fwaves1; spikes.waveforms_ch2=fwaves2;
        spikes.waveforms = [spikes.waveforms_ch1 spikes.waveforms_ch2];
        spikes.spiketimes=ftimes/1000;
        spikes.swtimes=fstimes/1000;

        spikes.fstimes=fstimes;
        spikes.ftimes=ftimes;

        spikes.sweep = datmat(:,2);
        spikes.trial = datmat(:,3);
        spikes.stimulus = stimcode(lp_f) * ones(size(spikes.ftimes));
        spikes.igorstim = datmat(:,5);

        spikes.Fs=32000;
        spikes.threshT=spikestarts;
        spikes.threshV=[-Inf, thrs];
        spikes.stimonset=stimonset;   %ms
        spikes.window=bckwindow;       %ms
        spikes.sweepd=sweepd;     %ms
        spikes.nsweeps = length(cuts)-1;

        newsweep = spikes.sweep + addsweep; newtrial = spikes.trial + addtrial;
        fstimes=(10000*newsweep) + ftimes;
        
%         if autosweepstim==1, sweepstim{lp_f}=newsweep; end
%         if autosweepall==1, sweepall=[sweepall; newsweep]; end
        
        if autosweepstim==1, sweepstim{lp_f}=Allsweepidx+addsweep; end
        if autosweepall==1, sweepall=[sweepall; Allsweepidx+addsweep]; end
        
        allspikes.waveforms_ch1=[allspikes.waveforms_ch1;fwaves1];       
        allspikes.waveforms_ch2=[allspikes.waveforms_ch2;fwaves2];

        allspikes.spiketimes=[allspikes.spiketimes;ftimes/1000];
        allspikes.swtimes=[allspikes.swtimes;fstimes/1000];
        allspikes.fstimes=[allspikes.fstimes;fstimes];
        allspikes.ftimes=[allspikes.ftimes;ftimes];
        allspikes.sweep = [allspikes.sweep; newsweep];
        allspikes.trial = [allspikes.trial; newtrial];


        allspikes.stimulus = [allspikes.stimulus; stimcode(lp_f) * ones(size(spikes.ftimes))];
        allspikes.igorstim = [allspikes.igorstim; datmat(:,5)];
        allspikes.nsweeps = allspikes.nsweeps + spikes.nsweeps;
        if consweep==1, addsweep = addsweep + spikes.nsweeps; addtrial = addtrial + max(spikes.trial); end

        save ([ name1 '_stereothr' [num2str(round(sdcuts(lpthrs)))]], ...
            'thrs', 'sweepd', 'spikes', 'stimonset','ndatmat','spikestarts')

        clear tempdat temp temp2 fwaves1 fwaves2 tempthr g_sp g_spt spikes

   %     close all; 
        lpthrs
    end % threshold loop

end % file loop

save temp
if (exist('sweepall')), allspikes.sweepall=sweepall; end
if (exist('sweepstim'))
    for i=1:length(unqstim)
        cmd=sprintf('allspikes.sweep_%d = sweepstim{i};',unqstim(i)); eval(cmd);
    end
end

allspikes.waveforms=[allspikes.waveforms_ch1 allspikes.waveforms_ch2];
allspikes.Alligorstim=Alligorstim;
spikes = allspikes;
spikes.Fs=32000;
spikes.threshT=spikestarts;
spikes.threshV=[-Inf, thrs];
spikes.stimonset=stimonset;   %ms
spikes.window=bckwindow;       %ms
spikes.sweepd=sweepd;     %ms

aspkfile = (['all_' recname '_thr' [num2str(round(sdcuts(lpthrs)))]]);

save ([ recname 'stereo_all_thr' [num2str(round(sdcuts(lpthrs)))]], 'spikes');