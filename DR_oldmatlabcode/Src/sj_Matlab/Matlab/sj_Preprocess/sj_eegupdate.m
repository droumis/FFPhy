function sj_eegupdate(prefix,days, epochs, idealtet, updatetets)
% Shantanu - Dec 2012. For EEGs o new DSP, usually eeg number >16, need to align to starttime
% of other eegs and trim it

% sj_eegupdate('HPb',[6], [1:5],15,16:20);

% sj_eegupdate('HPb',[3], [1:5],15,17:20);
% sj_eegupdate('HPb',[2], [1:5],7,17:20);
% sj_eegupdate('HPa',[3:4], [1:5],1,17:20);
% sj_eegupdate('HPa',[2], [1:5],1,17:20);
% sj_eegupdate('HPb',[1], [1:7],7,15:20);


if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epochs = 1:5; %
end
if nargin<4,
    idealtet = 1; %
end
if nargin<5,
    updatetets = 17:21; %
end


% SET DATA
% -------------------------------------------

switch prefix
    case 'HPa'
        rawdir = '/data25/sjadhav/HPExpt/HPa/';
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
    case 'HPb'
        rawdir = '/data25/sjadhav/HPExpt/HPb/';
        directoryname = '/data25/sjadhav/HPExpt/HPb_direct/';
    case 'HPc'
        rawdir = '/data25/sjadhav/HPExpt/HPc/';
        directoryname = '/data25/sjadhav/HPExpt/HPc_direct/';
end
dir2=directoryname;
savecopydir = [directoryname,'EEG/OriEEG/'];


if (idealtet<10)
    tetstring = ['0',num2str(idealtet)];
else
    tetstring = num2str(idealtet);
end

for d = 1:length(days)
    day = days(d);
    if (day<10)
        daystring = ['0',num2str(day)];
    else
        daystring = num2str(day);
    end
    
    for e = 1:length(epochs),
        
        epoch = epochs(e);
        % Get EEG for current epoch for the ideal tet
        % -------------------------------------------
        cd([directoryname,'/EEG/']);
        curreegfile = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',tetstring];
        load(curreegfile);
        %lfp = eeg{day}{epoch}{tet}.data;
        starttime = eeg{day}{epoch}{idealtet}.starttime; % This is in secs
        %endtime = starttime + (length(lfp)-1) * (1 / params.Fs);
        idealtimes = geteegtimes(eeg{day}{epoch}{idealtet});
        clear eeg
        
        for t = 1:length(updatetets)
            currtet = updatetets(t);
            if (currtet<10)
                ctetstring = ['0',num2str(currtet)];
            else
                ctetstring = num2str(currtet);
            end
            
            eegfile = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',ctetstring];
            load(eegfile);
            currtimes = geteegtimes(eeg{day}{epoch}{currtet});
            
            % Proceed only if lengths are actually mis-matched. A check.
            
            if length(currtimes)~=length(idealtimes);
               disp(['Length of curr tet is ',num2str(length(currtimes))]);
               disp(['Length of ideal tet is ',num2str(length(idealtimes))]);
                
                % Save a copy of this original in "OriEEG"
                savefile = [savecopydir,prefix,'eeg',daystring,'-',num2str(epoch),'-',ctetstring];
                save(savefile,'eeg');
                
                % Old Method - Was Assuming the bad tet had earlier starttimes and ok or later end times than ideal tet.
                % Have to check for all possibilities
                %                 stx = lookup(idealtimes(1),currtimes);
                %                 endx = lookup(idealtimes(end),currtimes);
                %                 enddiff = length(currtimes)-endx;
                %                 disp(['Tet',num2str(currtet),'Ep',num2str(e),':',' Start Diff in samples is ',num2str(stx)]);
                %                 disp(['Tet',num2str(currtet),'Ep',num2str(e),':',' End Diff in samples is ',num2str(enddiff)]);
                %                 %Trim EEG
                %                 eeg{day}{epoch}{currtet}.data = eeg{day}{epoch}{currtet}.data(stx:endx);
                %                 %And update starttime
                %                 eeg{day}{epoch}{currtet}.starttime = starttime;
                
                % Have to check for all possibilities
                
                stx = 1; endx = length(currtimes); % Default - No change
                stadd = 0; endadd = 0;
                
                
                % a) Check for trim or add at start
                
                if currtimes(1)~=idealtimes(1);
                    if currtimes(1) < idealtimes(1), % EEG starts earlier than ideal, so trim extra samples at start
                        stx = lookup(idealtimes(1),currtimes);
                        if stx==1, stx=stx+1;end
                        disp(['Tet',num2str(currtet),'Ep',num2str(e),' - Trim samples at Start: ',num2str(stx-1)]);
                    else
                        stadd = lookup(currtimes(1),idealtimes);
                        disp(['Tet',num2str(currtet),'Ep',num2str(e),' - Add samples at Start: ',num2str(stadd-1)]);
                    end
                end
                
                % b) Check for trim or add at end
                
                if currtimes(end) ~= idealtimes(end);
                    if currtimes(end) > idealtimes(end) % EEG ends later than ideal, so trim extra samples at end
                        endx = lookup(idealtimes(end),currtimes);
                        enddiff = length(currtimes)-endx;
                        disp(['Tet',num2str(currtet),'Ep',num2str(e),' - Trim samples at End: ',num2str(enddiff-1)]);
                    else
                        endid = lookup(currtimes(end),idealtimes);
                        endadd = length(idealtimes) - endid;
                        disp(['Tet',num2str(currtet),'Ep',num2str(e),' - Add samples at End: ',num2str(endadd)]);
                    end
                end
                
                % Check if TRIMMING is needed
                if stx~=1 || endx~=length(currtimes);
                    % Trim EEG
                    eeg{day}{epoch}{currtet}.data = eeg{day}{epoch}{currtet}.data(stx:endx);
                end
                
                % Check if lengths are equal. If not, check for add samples, or else its a sampling rate error
                if length(eeg{day}{epoch}{currtet}.data) ~= length(idealtimes),
                    disp(['Lengths still not equal after trim']);
                    
                    if stadd~=0 || endadd~=0
                        disp(['Its ok. Need to add samples']);
                        
                        if stadd~=0
                            staddvec = eeg{day}{epoch}{currtet}.data(1)*ones(stadd-1,1);
                        else
                            staddvec=[];
                        end
                        
                        if endadd~=0
                            endaddvec = eeg{day}{epoch}{currtet}.data(end)*ones(endadd,1);
                        else
                            endaddvec=[];
                        end
                        
                        eeg{day}{epoch}{currtet}.data = [staddvec;eeg{day}{epoch}{currtet}.data;endaddvec];
                        
                    end   % End add samples if
                    
                end % End check length if
                
                % Check again after adding samples
                if length(eeg{day}{epoch}{currtet}.data) ~= length(idealtimes),
                    disp(['Error - Lengths still do not match for Day',num2str(day), ' Ep',num2str(epoch),' Tet', num2str(currtet)]);
                    disp(['So its a Sampling rate error - needs interpolation. To be implemented']);
                    keyboard;
                end
                
                % all is ok. Update starttime to ideal starttime and save in EEG directory. this will overwrite existing eeg file
                eeg{day}{epoch}{currtet}.starttime = starttime;
                disp(['All ok after changes. Saving updated eegfile: ',eegfile]);
                disp(['Length of curr tet is ',num2str(length(eeg{day}{epoch}{currtet}.data))]);
                disp(['Length of ideal tet is ',num2str(length(idealtimes))]);
                save(eegfile,'eeg');
                
            else
                
                disp(['Length is fine for Tet',num2str(currtet),'Ep',num2str(epoch),'. No changes implemented']);
                
            end
            
        end  % end tets
        
        
    end % end epochs
    
    
end % end days






