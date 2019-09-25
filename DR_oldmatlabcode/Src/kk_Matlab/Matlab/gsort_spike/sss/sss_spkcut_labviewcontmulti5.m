function [spikes,aspkfile]=sss_spkcut_labviewcontmulti5(files, ch, Fs, sdcuts, threshold, stimcode, params, consweep, sweepall, sweepstim)
%sss_spkcut_labviewcontmulti5('spike.dat', [1:4], 32000, 4, 1);
%VER 5: JUN 23, 2006
%CHANGE CUMULATIVE SPIKE COUNT

%% threshold=0, +ve spikes; threshold=1, -ve spikes

%%%% LOAD Labview, OR USE DATA&TIME FIELDS FROM DAQREAD, AND CUT SPIKES AND MAKE A SPIKES STRUCTURE
%%%% KEEPING ALL FIELDS FOR BACKWARD COMPATIBILITY: NOT NECESSARY. THIS IS
%%%% NOT SWEEP BASED ACQUISITION
% sss_spkcut_labviewcontmulti5({'spike.dat'},[1:4],32000,4,1);
% sss_spkcut_labviewcontmulti5({'spike.dat'},[5:8],32000,4,1);
% sss_spkcut_labviewcont4({'x1';'x2'},[1:4],30000,5,0,0)
% sss_spkcut_labviewcont4({'st1e1_r1';'st1e1_r2'}, [1:4],30000,4,0,0,
% params, [0 1 0 1 100], 1, [1:1090], {[1:200 601:800]; [201:600 801:1000]; [1001:1090]});

%%% READ THIS SECTION "CAREFULLY!" FOR FIRST TIME USERS

%%% Modified: Shantanu Jadhav, Nov21, 2005

% This is a function to be used as a callback in GSORTM,
% and is written to automatically detect and cut spikes at various thresholds.
% Threshold calculation is in respect to the background noise level before the
% stimulus onset.  The user specifies one standard deviation value and the program
% calculates the mean voltage values and standard deviations.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MANDATORY %%%%%%%%%%%%%%%%%%%%%%%%%%%

%   FILES is a cell array with the .swp file names written individually in seperate
%   lines. ex. {'r49g5P1'; 'r49g5R2';'r49g5SR1'}.  PLEASE NOTE THAT THE FIRST FILE
%   SHOULD BE THE ONE WHERE THE THRESHOLD INFORMATION WILL BE CALCULATED.

%%%%%%%%%%%%%% OPTIONAL, IF NO INPUT, DEFAULTS WILL BE USED:%%%%%%%%%%%%%%
%%%% FIELDS necessary for SSS_STIM_ARCHIVE are PARAMS, STIMCODE(SPIKES.STIMULUS),
%%%% SWEEPALL(SPIKES.SWEEPALL), SWEEPSTIM(SPIKES.SWEEPSTIM) %%%%%%%%%%%%%%

%   FILES: Labview Files with continuous data: Cut into 5 min segments and loaded by
%   labview_loadspike4: Spikes cut from each segment and then spliced
%   together in the "*all*" file
%
%   CH: Channel Numbers in Labview DATA
%
%   Fs: Sampling Frequency
%
%   SDCUTS is the values to be used to calculate the thresholds.  For example if as
%   SDCUTS there is only one value, say 3, then the program will calculate the mean
%   plus the three sds as the threshold. Example entries: 3|[3:5]|[3 5 9].
%   USE ONLY ONE VALUE
%
%   THRESHOLD = 0:+ve spikes; =1, -ve spikes;
%
%   STIMCODE: %% Ignore for continuous files. This is for backward compatibility. Use 0 
%  {(Tricky) For a super-structure stimulus code, which is higher in
%   hierarchy than igorstim. You need one stimcode for each .swp file input. The
%   default stimcodes are 0,1,2...length(files)-1. If you need multiple
%   stimcodes in one file, NEED TO MANUALLY PUT SWEEPALL AND SWEEPSTIM
%   In future, will have an IGORFIELD for this FOR EACH TRIAL, AND THAT WILL BE
%   DEFAULT. THERE WILL HAVE TO CHANGES TO DEFAULT CALCS OF SWEEPALL AND SWEEPSTIM
%   THEN.}
%
%   PARAMS: IGNORE- Backward Compatibility
%        stimonset, bckwindow, sweepd
%       STIMONSET is the onset of the whisker deflection in ms resolution.
%       BCKWINDOW is the window size where the average voltage and the standard
%       deviation of the voltage will be calculated. Resolution in ms.
%       SWEEPD is the sweep duration in ms resolution.
%       params can be empty

%   CONSWEEP: IGNORE- Backward Compatibility 
%   Default=1 implies sweeps will be pushed for multiple files (almost
%   always the case). You want a unique sweep number for all sweeps (for sorting
%   and isi reasons), each with a igorstim (which is in the SPIKES.ALLIGORSTIM field)
%
%   SWEEPALL: IGNORE- Backward Compatibility
%   Range of sweeps, 1st to last across all files. Default should
%   be done auto by adding sweeps from files by concatanating them
%
%   SWEEPSTIM: IGNORE- Backward Compatibility
%   Sweep range for each stimulus: one file => one stimulus =>
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
format long

%% SHIFT CHANNELS %%%%%%%%%%
%ch=ch+1;
%% SHIFT CHANNELS %%%%%%%%%%

if (nargin < 3), sfreq=32; end % sampling frequency in kHz.
if (nargin < 4), sdcuts=4; end
if (nargin < 5), threshold=0; end
if (nargin < 6), stimcode=(1:size(files,1))-1; end %% default stimcode = 0,1,2...lth(files)-1
if (nargin < 7),
    stimonset=200; bckwindow=100; sweepd=1000;
else
    if (isempty(params.stimonset)), stimonset = 200; else, stimonset=params.stimonset; end
    if (isempty(params.bckwindow)), bckwindow = 100; else, bckwindow=params.bckwindow; end
    if (isempty(params.sweepd)), sweepd = 1000; else, sweepd=params.sweepd; end
end
if (nargin < 8), consweep=1; end

if (exist('sweepall')),
    if (isempty(sweepall)),autosweepall=1; sweepall=[];else,autosweepall=0;end
else
    autosweepall=1; sweepall=[];
end

if (exist('sweepstim')),
    if (isempty(sweepstim)),autosweepstim=1; sweepstim=[];else,autosweepstim=0;end
else
    autosweepstim=1; sweepstim=[];
end
%if (exist('sdcuts')==0), sdcuts=5; end
if length(stimcode) < size(files,1)
    stimcode(end+1:size(files,1))=0;
elseif length(stimcode) > size(files,1)
    stimcode(size(files1,1):end)=[];
end

 %%%%%%%%%%%%%%%%
    % Create variables to use later in the code
    longsp=0; bsize=1;
    if (threshold==0), spikestarts=14; else, spikestarts=10; end
 %%%%%%%%%%%%%%%%%


unqstim=unique(stimcode);
recname = 'spike';
addsweep=0; addtrial=0; maxsweep=0;

allspikes.ftimes=[]; allspikes.fstimes=[];

stop=0;
lp_f=0;
while stop==0
    lp_f=lp_f+1;
    %Get data and time from channels, can also get abstime

    name1 = 'spike';
    fn='spike.dat';

    %%% DIVIDE INTO 5MIN SEGMENTS FOR LOADING;  
    t1=300*(lp_f-1); t2=300*(lp_f);
    lp_f
    t1, t2
    %data=labview_loadspike4(fn,length(ch),32000,t1,t2);
    
    %nch=length(ch);
    nch=4;
    
    data=labview_loadspike4(fn,nch,Fs,t1,t2);  %% NEED TO READ ALL 4 channels
    data=data(ch,:); data=data';
   
    
    %%%%% REMOVE HIGH NOISE EPOCHS: CHECK FEASIBILITY %%%%%%
    %    
            x=find((data(:,2)<-0.48)|(data(:,2)>0.48));
            x(find(x<=51))=[]; x(find(x>=length(data)-51))=[]; 
            for i=1:length(x), data(x(i)-50:x(i)+50,:)=0;end
    %    
     %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
     
    %% THRESHOLD CROSSING ON FIRST CHANNEL ONLY: AND for first file only
    if lp_f==1 % if the first cycle then calculate the threshold
        tempthr= [mean(mean(data(1:Fs,2))) mean(std(data(1:Fs,2)))];
        for lpsdc=1:length(sdcuts)
            if (threshold==0),
                thrs(lpsdc)= tempthr(1)+sdcuts(lpsdc)*tempthr(2);
            else
                thrs(lpsdc)= tempthr(1)-sdcuts(lpsdc)*tempthr(2);
            end
            thrs
        end
    end

    thrs=-0.08

    for lpthrs=1:length(thrs);
        thr= thrs(lpthrs);      
        for nch=1:length(ch)
            %cmd=sprintf('spikes.waveforms_ch%d = [];',ch(nch)); eval(cmd);
            cmd=sprintf('spikes.waveforms_ch%d = [];',nch); eval(cmd);
        end
        spikes.fstimes=[]; spikes.ftimes=[];
        h = waitbar(0,'Making spiketimes');

        %% SAVE IGORSTIM NUMBER IRRESPECTIVE OF SPIKE OR NOT %%%
        switch length(ch),

            case 2,
                indices = find( (data(:,1)<=thr) | (data(:,2)<=thr)); % values over threshold
            case 3,
                indices = find( (data(:,1)<=thr) | (data(:,2)<=thr) | (data(:,3)<=thr));
            case 4,
                indices = find( (data(:,1)<=thr) | (data(:,2)<=thr) | (data(:,3)<=thr) | (data(:,4)<=thr));
            otherwise
                indices = find( (data(:,1)<=thr) | (data(:,2)<=thr) | (data(:,3)<=thr) | (data(:,4)<=thr));
        end

        thr_index = indices([find(diff(indices)>1)+1]); % threshold crossings in +/- direction: point of threshold cross
        thr_index =[indices(1); thr_index]; %spikes start:  index of threshold crossing

        %thrdiff=find(diff(thr_index)<=8); %% Corresponds to 0.25ms diff; discard spikes within this range? Multiple counting?
        thrdiff=find(diff(thr_index)<=5);
        thr_index(thrdiff+1)=[];

        thr_index=thr_index(find ((thr_index>spikestarts+1) & (thr_index<length(data)-(32-spikestarts)))); %Refine index of
        % threshold crossing by getting those spikes that are entirely within the 1 ms: 32 pt window

        spktime= round(thr_index * (10^4/Fs)); %  constant converts the spiketimes from Fs sampling to  resolution of 0.1 ms       
        spktime= spktime/10; %CONVERT TO ms format: SHOULD ROUND IT TO TO JUST ONE DECIMAL AFTER .

        %%MATLAB SPKTIME
        %             spktime = time(thr_index); %Time of spike: thrs crossing: resolution of 1/Fs, in sec
        %             spktime=spktime*(10^3); %%CONVERT TO ms: SHOULD ROUND IT TO TO JUST ONE DECIMAL AFTER .
        %             spktime=roundn(spktime,-1); %% Resolution of 0.1 ms

        %spktime = round(spktime *(10^4)/Fs);    % constant converts the spiketimes from Fs sampling to 0.1 ms

        % 1 ms: or 32 pt spikes
        for sp=1:length(thr_index);
            spikes.ftimes=[spikes.ftimes; spktime(sp)];
            waitbar(sp/length(thr_index));
        end       
        
        %%%% CHANGED ON MAR28, 2006: USE THIS AS A CUMULATIVE COUNT OF SPIKETIME IN MSEC
        %%%% & SECONDS RESPECTIVELY%%%%%%%
        %spikes.fstimes=spikes.ftimes + 1000*(addsweep+1);     
        spikes.fstimes=spikes.ftimes + 1000*(addsweep);    %% DO NOT PUSH BY 1000MS: JUN 23. 2006
        
        %%%%%%%%%% CHANGED ON MAR28, 2006 %%%%%%%%%%%%%     
        %spikes.sweep=floor((spikes.ftimes)/1000)+1+addsweep;    %%% DIVIDE INTO 1 SEC ARTIFICIAL SWEEPS
        
        spikes.sweep=floor((spikes.ftimes)/1000)+addsweep;        
        maxsweep=max(spikes.sweep);
        
        
        
        spikes.Fs=Fs; spikes.threshT=spikestarts; spikes.threshV=[-Inf, thrs];
        
        allspikes.ftimes=[allspikes.ftimes;spikes.ftimes];
        allspikes.fstimes=[allspikes.fstimes;spikes.fstimes];  %%% UPDATED FSTIMES
        
        if consweep==1, 
            if size(data,1)<Fs*300      %%% LAST FILE ADDSWEEP WILL GIVE THE SECOND FOR LAST SPIKE TIME 
                addsweep = maxsweep+1;
            else
                addsweep = addsweep + 300;
            end                    
        end
        save ([ name1 '_multi_' [num2str(lp_f)] '_tet_thr' [num2str(thrs*1000)]], 'spikes')
        clear fstimes indices thr_index spktime spikes time
        
    end % thrs loop
    
    addsweep
    if size(data,1)<Fs*300
        stop=1;
    end   
    clear data
    
end % file loop

save temp
spikes=allspikes;

spikes.Fs=Fs; spikes.threshT=spikestarts; 
if thrs<0, spikes.threshV=[thrs Inf]; else, spikes.threshV=[-Inf thrs];end
spikes.nsweeps=floor(max(spikes.fstimes)/1000);

spikes.hierarchy.assigns=ones(size(spikes.fstimes));

save ([ recname 'multi_all_thr' [num2str(thrs*1000)]], 'spikes');
