function [spikes,aspkfile]=sss_spkcut_matlabcont2(files, ch, Fs, sdcuts, stimcode, params, consweep, sweepall, sweepstim)

%%%% LOAD MATLAB DAQ FILE, OR USE DATA&TIME FIELDS FROM DAQREAD, AND CUT SPIKES AND MAKE A SPIKES STRUCTURE
%%%% KEEPING ALL FIELDS FOR BACKWARD COMPATIBILITY: NOT NECESSARY. THIS IS
%%%% NOT SWEEP BASED ACQUISITION
% sss_spkcut_matlabcont1({'st1e1_r1';'st1e1_r2}', [4:7],30000,4)
% params, [0 1 0 1 100], 1, [1:1090], {[1:200 601:800]; [201:600 801:1000]; [1001:1090]});

%%% READ THIS SECTION "CAREFULLY!" FOR FIRST TIME USERS

%%% Modified in Feb2006 from sss_autospikecut (Shantanu Jadhav, Nov21,
%%% 2005)

% SSS_AUTOSPIKECUT_1MS_BASIC is a function to be used as a callback in GSORTM,
% and is written to automatically detect and cut spikes at various thresholds.
% Threshold calculation is in respect to the background noise level before the
% stimulus onset.  The user specifies one standard deviation value and the program
% calculates the mean voltage values and standard deviations.

% THIS FILE CUTS 1 MS SPIKES FOR SINGLE CHANNEL ONLY (32 samples at Fs kHz)
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
format long
% Create variables to use later in the code
longsp=0; bsize=1; spikestarts=14;
%% SHIFT CHANNELS %%%%%%%%%%
ch=ch+1;
%% SHIFT CHANNELS %%%%%%%%%%

if (nargin < 3), sfreq=30; end % sampling frequency in kHz.
if (nargin < 4), sdcuts=4; end
if (nargin < 5), stimcode=(1:size(files,1))-1; end %% default stimcode = 0,1,2...lth(files)-1
if (nargin < 6),
    stimonset=200; bckwindow=100; sweepd=1000;
else
    if (isempty(params.stimonset)), stimonset = 200; else, stimonset=params.stimonset; end
    if (isempty(params.bckwindow)), bckwindow = 100; else, bckwindow=params.bckwindow; end
    if (isempty(params.sweepd)), sweepd = 1000; else, sweepd=params.sweepd; end
end
if (nargin < 7), consweep=1; end

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
if length(stimcode) < size(files,1)
    stimcode(end+1:size(files,1))=0;
elseif length(stimcode) > size(files,1)
    stimcode(size(files1,1):end)=[];
end

unqstim=unique(stimcode);
recname = strtok(files{1},'.');
% [base1,base2] = strtok(files{1},'L');
% recname=[base1 base2(1:end-1)]; %    Rec Site Name

addsweep=0; addtrial=0;
allspikes.waveforms=[];
for nch=1:length(ch)
    %cmd=sprintf('allspikes.waveforms_ch%d = [];',ch(nch)); eval(cmd);
    cmd=sprintf('allspikes.waveforms_ch%d = [];',nch); eval(cmd);
end
%allspikes.waveforms_ch1=[];allspikes.waveforms_ch2=[];
%allspikes.waveforms_ch3=[];allspikes.waveforms_ch4=[];
allspikes.spiketimes=[]; allspikes.swtimes=[];
allspikes.fstimes=[];allspikes.ftimes=[];
allspikes.sweep = [];allspikes.trial = [];
allspikes.stimulus = []; allspikes.igorstim = [];allspikes.nsweeps = 0;
Alligorstim=[]; Allsweepidx=[];

for lp_f=1:size(files,1)

    %Get data and time from channels, can also get abstime
    name1 = strtok(files{lp_f},'.');
    [data time]=daqread([files{lp_f}],'Channels',ch);

    %% THRESHOLD CROSSING ON FIRST CHANNEL ONLY: AND for first file only
    if lp_f==1 % if the first cycle then calculate the threshold
        %tempdat=data(:,8+(stimonset*sfreq-bckwindow*sfreq):8+(stimonset*sfreq));
        tempthr= [mean(mean(data(1:30000,1))) mean(std(data(1:30000,1)))];
        for lpsdc=1:length(sdcuts)
            thrs(lpsdc)= tempthr(1)+sdcuts(lpsdc)*tempthr(2);
        end
    end
    
    %thr=40;
    
    for lpthrs=1:length(thrs);
        thr= thrs(lpthrs);
        shape=[];
        for nch=1:length(ch)
            %cmd=sprintf('spikes.waveforms_ch%d = [];',ch(nch)); eval(cmd);
            cmd=sprintf('spikes.waveforms_ch%d = [];',nch); eval(cmd);
        end
        spikes.fstimes=[]; spikes.ftimes=[];
        h = waitbar(0,'Creating the waveforms');

        %% SAVE IGORSTIM NUMBER IRRESPECTIVE OF SPIKE OR NOT %%%
        Alligorstim = [Alligorstim lp_f];
        Allsweepidx = [Allsweepidx lp_f];  %%% Save sweep no irrespective of spike or not in sw For sweepall & sweepsti
        %% Keeping these fields for backward compatibility %%

        % Go and find the threshold passing events on Ch1 of array ch: Thrs on this 1st channel only
        
        %%%%%%THRESHOLDING ON ALL 4 CHANNELS  %%%%%%%%%%%%%%%%%%%%%%%
        
        %if length(find(data(:,1)>=thr))>=1
            
              %indices = find (data(:,1)<=thr); % values over threshold
            indices = find( (data(:,1)>=thr) | (data(:,2)>=thr) | (data(:,3)>=thr)|(data(:,4)>=thr) ); % values over threshold
            thr_index = indices([find(diff(indices)>1)+1]); % threshold crossings in +/- direction: point of threshold cross
            thr_index =[indices(1); thr_index]; %spikes start:  index of threshold crossing

%             thrdiff=find(diff(thr_index)<=8); %% Corresponds to 0.25ms diff; discard spikes within this range? Multiple counting?
            thrdiff=find(diff(thr_index)<=5);
            thr_index(thrdiff+1)=[];
            
            thr_index=thr_index(find ((thr_index>spikestarts+1) & (thr_index<length(data)-(32-spikestarts)))); %Refine index of
            % threshold crossing by getting those spikes that are entirely within the 1 ms: 32 pt window
                       
            spktime = time(thr_index); %Time of spike: thrs crossing: resolution of 1/Fs, in sec

            spktime=spktime*(10^3); %%CONVERT TO ms: SHOULD ROUND IT TO TO JUST ONE DECIMAL AFTER .
            spktime=roundn(spktime,-1); %% Resolution of 0.1 ms
            %spktime = round(spktime *(10^4)/Fs);    % constant converts the spiketimes from Fs sampling to 0.1 ms

            % 1 ms: or 32 pt spikes
            for sp=1:length(thr_index);
                spikes.ftimes=[spikes.ftimes; spktime(sp)]; 
                %fspiketime=(lp_f*100000)+spktime(sp); %spikes.fstimes=[spikes.fstimes; fspiketime];
                for nch=1:length(ch)
                    %shape= data( thr_index(sp)-spikestarts:thr_index(sp)+(31-spikestarts),ch(nch) );
                    shape= data( thr_index(sp)-spikestarts:thr_index(sp)+(31-spikestarts),nch );
                    %cmd=sprintf('spikes.waveforms_ch%d = [spikes.waveforms_ch%d shape];',ch(nch), ch(nch)); eval(cmd);
                    cmd=sprintf('spikes.waveforms_ch%d = [spikes.waveforms_ch%d shape];',nch, nch); eval(cmd); clear shape
                end
                waitbar(sp/length(thr_index));
            end
        %end %end if


        name = strtok(files{1},'.'); spikes.waveforms=[];
        for nch=1:length(ch),
            %cmd=sprintf('spikes.waveforms = [ spikes.waveforms spikes.waveforms_ch%d];',ch(nch)); eval(cmd);
            cmd=sprintf('spikes.waveforms = [ spikes.waveforms; spikes.waveforms_ch%d];',nch); eval(cmd);
        end
        %spikes.waveforms = spikes.waveforms';
        
        %%%NOTE THAT HERE I AM ASSUMING THAT I WILL GET THE NECESSARY MAX
        %%%EXPONENT (FOR PUSHING TO FSTIMES BASED ON FILE NUMBER) FROM THE 1ST FILE ITSELF
        %%%THIS IS OK AS LONG AS FURTHER FILES DO NOT HAVE
        %%%ORDER-PF-MAGNITUDE-HIGHER SPIKETIMES
        if lp_f==1
            for i=1:10, x=max(spikes.ftimes)/(10^i); if x<1, exp=i; break; end; end
        end
        expo(lp_f)=exp;
        spikes.fstimes=(lp_f)*(10^exp)+spikes.ftimes;
        spikes.spiketimes=spikes.ftimes/1000; spikes.swtimes=spikes.fstimes/1000;

        spikes.sweep = lp_f*ones(size(spikes.ftimes)); spikes.trial = lp_f*ones(size(spikes.ftimes));
        spikes.stimulus = stimcode(lp_f) * ones(size(spikes.ftimes)); spikes.igorstim = stimcode(lp_f) * ones(size(spikes.ftimes));

        spikes.Fs=Fs; spikes.threshT=spikestarts; spikes.threshV=[-Inf, thrs];
        spikes.stimonset=stimonset; spikes.window=bckwindow; spikes.sweepd=sweepd;
        spikes.nsweeps = size(files,1);

        newsweep = spikes.sweep + addsweep; newtrial = spikes.trial + addtrial;
        fstimes=(10000*newsweep) + spikes.ftimes;
        if autosweepstim==1, sweepstim{lp_f}=Allsweepidx+addsweep; end
        if autosweepall==1, sweepall=[sweepall Allsweepidx+addsweep]; end

        for nch=1:length(ch),
            %cmd=sprintf('allspikes.waveforms_ch%d = [allspikes.waveforms_ch%d  spikes.waveforms_ch%d];',ch(nch), ch(nch)); eval(cmd);
            cmd=sprintf('allspikes.waveforms_ch%d = [allspikes.waveforms_ch%d  spikes.waveforms_ch%d];',nch, nch, nch); eval(cmd);
        end

        allspikes.spiketimes=[allspikes.spiketimes;spikes.ftimes/1000]; allspikes.ftimes=[allspikes.ftimes;spikes.ftimes];
        allspikes.swtimes=[allspikes.swtimes;spikes.fstimes/1000]; allspikes.fstimes=[allspikes.fstimes;spikes.fstimes];  %%% UPDATED FSTIMES
        allspikes.sweep = [allspikes.sweep; newsweep]; allspikes.trial = [allspikes.trial; newtrial];
        allspikes.stimulus = [allspikes.stimulus; stimcode(lp_f) * ones(size(spikes.ftimes))];
        allspikes.igorstim = [allspikes.igorstim; stimcode(lp_f) * ones(size(spikes.ftimes))];
        allspikes.nsweeps =  size(files,1);
        if consweep==1, addsweep = addsweep + 1; addtrial = addtrial + 1; end

        save ([ name1 'cut_tet_thr' [num2str(round(sdcuts(lpthrs)))]], 'spikes')
        clear shape fstimes indices thr_index spktime spikes data time
        %lpthrs
    end % threshold loop
    lp_f

end % file loop

save temp
if (exist('sweepall')), allspikes.sweepall=sweepall; end
if (exist('sweepstim'))
    for i=1:length(unqstim), cmd=sprintf('allspikes.sweep_%d = sweepstim{i};',unqstim(i)); eval(cmd); end
end
allspikes.Alligorstim=Alligorstim;
for nch=1:length(ch),
    %cmd=sprintf('allspikes.waveforms = [ allspikes.waveforms allspikes.waveforms_ch%d];',ch(nch)); eval(cmd);
    cmd=sprintf('allspikes.waveforms = [ allspikes.waveforms; allspikes.waveforms_ch%d];',nch); eval(cmd);
end
allspikes.waveforms = allspikes.waveforms';
%for nch=1:length(ch), cmd=sprintf('allspikes.waveforms_ch%d = [allspikes.waveforms_ch%d];'nch, nch); eval(cmd); end

spikes=allspikes;

%for i=1:10, x=max(spikes.ftimes)/(10^i); if x<1, exp=i; break; end; end
%exp=max(expo);
%spikes.fstimes=(10^exp)+spikes.ftimes;
spikes.Fs=Fs; spikes.threshT=spikestarts; spikes.threshV=[-Inf, thrs];
spikes.stimonset=stimonset; spikes.window=bckwindow;spikes.sweepd=sweepd;
spikes.nsweeps=floor(max(spikes.ftimes)/spikes.sweepd);

aspkfile = (['all_' recname '_thr' [num2str(round(sdcuts(lpthrs)))]]);
save ([ recname '_cut_tet_all_thr' [num2str(round(sdcuts(lpthrs)))]], 'spikes');

siz=size(spikes.waveforms);
spikes.waveforms=reshape(spikes.waveforms,[siz(2) siz(1)]);
siz=size(spikes.waveforms_ch1);
for nch=1:length(ch),
    cmd=sprintf('spikes.waveforms_ch%d = reshape(spikes.waveforms_ch%d,[siz(2) siz(1)]);',nch,nch); eval(cmd);
end

save ([ recname '_resh_cut_tet_all_thr' [num2str(round(sdcuts(lpthrs)))]], 'spikes');

