function [spikes,aspkfile]=sss_autospikecut (files,stimcode,consweep)

% cr_autospikecut ({'CR10par3_e1L4_AllWh1';'CR10par3_e1L4_AllWh2'}, [0 0], 200, 100, 6, 700, 500)
% cr_autospikecut ({'A27p3r1_TC';'A27p3r1_1';'A27p3r1_2';'A27p3r1_3';'A27p3r1_4'}, [0 1 2 3 4], 200, 100, [11.5], 1000, 50, 390);
% Assuming same no of sweeps for all stimuli 


% AUTOSPIKECUT is a function in the NEWIGOR2MAT toolbox, and is written to automatically detect
% and cut spikes at various thresholds.  Threshold calculation is in respect to the background 
% noise level before the stimulus onset.  The user specifies one or more standard deviation values
% and the program calculates the mean voltage values and standard deviations. 
% 
% INPUTS:
%   FILES is a cell array with the file names written individually in seperate lines.  
%       ex. {'r49g5P1'; 'r49g5R2';'r49g5SR1'}.  PLEASE NOTE THAT THE FIRST FILE SHOULD BE THE 
%       ONE WHERE THE THRESHOLD INFORMATION WILL BE CALCULATED.
%   STIMONSET is the onset of the whisker deflection in ms resolution.
%   BCKWINDOW is the window size where the average voltage and the standard deviation of the
%       voltage will be calculated. Resolution in ms. 

%   SDCUTS is the values to be used to calculate the thresholds.  For example if as SDCUTS there
%       is only one value, say 3, then the program will calculate the mean plus the three standard
%       deviation as the threshold. Example entries: 3 | [3:5] | [3 5 9].   
%   SWEEPD is the sweep duration in ms resolution.
%   WSIZE is the duration of the background activity desired to be studied.  
%   N_OF_SWEEPS how many sweeps you have recorded for a given whisker's deflection (i.e. 900).

% my_autospikecut ({'r49g5P1'; 'r49g5R2';'r49g5SR1'}, 100, 50, [4:5], 600, 50, 100)

% Tansu CELIKEL
% All rights reserved

warning off;
stimonset = 200;
bckwindow=100;
sdcuts=5;
sweepd=700;
if length(stimcode) < size(files,1)
    stimcode(end+1:size(files,1))=0;
elseif length(stimcode) > size(files,1)
    stimcode(size(files,1):end)=[];
end

% Create variables to use later in the code
longsp=0; bsize=1; spikestarts=14;
sfreq=32; % sampling frequency in kHz. 
[base1,base2] = strtok(files{1},'L');
recname=[base1 base2(1:end-1)]; %    Rec Site Name

addsweep=0; addtrial=0;
allspikes.waveforms=[];
allspikes.spiketimes=[];
allspikes.swtimes=[];

allspikes.fstimes=[];
allspikes.ftimes=[];

allspikes.sweep = [];
allspikes.trial = [];
allspikes.stimulus = [];
allspikes.igorstim = [];
allspikes.nsweeps = 0;



for lp_f=1:size(files,1)
    
    name = ([files]);
    data=load ([files{lp_f}]);
    
    %data = inh;
    h = waitbar(0,'Reorganizing the data into matrix...');
    
    data = [-9999 data' -9999];
    
    cuts = find (data == -9999);
    datmat2 = zeros (length(cuts)-1,(cuts(2)-cuts(1))+1); % checked
    
    for lp_c = 1: (length (cuts)-1)
        datmat2 (lp_c,:)= data (cuts(lp_c):cuts(lp_c+1)); %  checked. for an example data set:see test.xls
        waitbar(lp_c/(length (cuts)-1));
        a=1;
    end
    
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
    
    close force; clear data ndata h lp_c
    
    % CAN CALCULATE THRESHOLD FOR ALL FILES SEPARATELY: KEEP FOR NOW
    
    sdcuts=6;
    
    if lp_f==1 % if the first cycle then calculate the threshold
        tempdat=datmat2(:,8+(stimonset*sfreq-bckwindow*sfreq):8+(stimonset*sfreq));
        tempthr= [mean(mean(tempdat)) mean(std(tempdat))];
        
        for lpsdc=1:length(sdcuts)
            thrs(lpsdc)= tempthr(1)+sdcuts(lpsdc)*tempthr(2);    
        end    
    end
    
    
    % thrs = 9;
    
    for lpthrs=1:length(thrs);
        
        thr= thrs(lpthrs);
        datmat=[]; fstimes=[]; fspikes=[]; indct=1;
        h = waitbar(0,'Creating the waveforms...');
        
        %keyboard
        % Time to find the waveforms and their corresponding spike times.
        % First create a for loop to go over the sweep data (sw).
        for lp_sw = 1: size (datmat2,1)
            
            %Create a temporary vector where you can write a single sweep data
            temp2 = datmat2(lp_sw,:);
            temp = temp2 (8:end-1); % clean the non A/D values in the temp
            % Go and find the threshold passing events
            if length(find(temp>=thr))>=1
                g_sp= find (temp>=thr); % values over threshold
                g2_sp= g_sp([find(diff(g_sp)>1)+1]); % threshold crossings in +/- direction
                g2_sp=[g_sp(1) g2_sp]; %spikes start
                               
                g2_sp=g2_sp(find ((g2_sp>spikestarts+1) & (g2_sp<size(temp2,2)-(40-spikestarts))));
                % get those spikes that are entirely within the 1 ms window
                
                g2_spt= round(g2_sp * 0.3125);    % 0.3125 constant converts the spiketimes from 32kHz sampling to 0.1 ms  
                
                %stc(lp_sw)=length(g2_spt);
                % 1 ms spikes
                for lp_sp=1:length(g2_sp);
                    temps= temp(g2_sp(lp_sp)-spikestarts:g2_sp(lp_sp)+(31-spikestarts));            
                    fspikes= [fspikes; temps];
                    fspiketime=(lp_sw*10000)+g2_spt(lp_sp); fstimes=[fstimes; fspiketime];  
                    datmat(indct,1)=-9999; % seperator 
                    datmat(indct,2)=temp2(1,4) + 1; % sweep number
                    datmat(indct,3)=temp2(1,5); % trial number
                    datmat(indct,4)=temp2(1,6); % repetition number
                    datmat(indct,5)=temp2(1,7); % stimulus number
                    datmat(indct,6)=g2_spt(lp_sp); % time stamp
                    indct=indct+1;
                end
            end %end if
            waitbar(lp_sw/size (datmat2,1));
        end    
        
        
        fwaves=fspikes; 
        datmat= [datmat fwaves datmat(:,1)]; 
        ftimes= [datmat(:,6)]; 
        
        scthreshold=sdcuts(lpthrs);
        fstudied= files{lp_f};
        
        clear fspikes 
        close force
        
        nsweeps = length (cuts)-1;
        ndatmat=datmat;
        name = strtok(files{1},'_');
        
        spikes.waveforms=fwaves;
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
        allspikes.waveforms=[allspikes.waveforms;fwaves];
        allspikes.spiketimes=[allspikes.spiketimes;ftimes/1000];
        allspikes.swtimes=[allspikes.swtimes;fstimes/1000];
        allspikes.fstimes=[allspikes.fstimes;fstimes];
        allspikes.ftimes=[allspikes.ftimes;ftimes];
        allspikes.sweep = [allspikes.sweep; newsweep];
        allspikes.trial = [allspikes.trial; newtrial];
        
        
        allspikes.stimulus = [allspikes.stimulus; stimcode(lp_f) * ones(size(spikes.ftimes))];
        allspikes.igorstim = [allspikes.igorstim; datmat(:,5)];
        allspikes.nsweeps = allspikes.nsweeps + spikes.nsweeps;
        if consweep==1
            addsweep = addsweep + spikes.nsweeps; addtrial = addtrial + max(spikes.trial);    
        end
        
        save ([ files{lp_f} '_thr' [num2str(round(sdcuts(lpthrs)))]], ...
            'fwaves', 'ftimes', 'thrs', 'fstimes', 'datmat', 'sweepd', 'spikes', 'stimonset','sweepd','ndatmat','spikestarts')  
        
        clear tempdat temp temp2 temps tempthr g_sp g_spt Ratiovals spikes
        
        close all; close force
        lpthrs
    end % threshold loop
    
end % file loop

spikes = allspikes;
spikes.Fs=32000;
spikes.threshT=spikestarts;
spikes.threshV=[-Inf, thrs];
spikes.stimonset=stimonset;   %ms
spikes.window=bckwindow;       %ms
spikes.sweepd=sweepd;     %ms

aspkfile = ([recname 'all_thr' [num2str(round(sdcuts(lpthrs)))]]);       

save ([ recname 'all_thr' [num2str(round(sdcuts(lpthrs)))]], 'spikes');