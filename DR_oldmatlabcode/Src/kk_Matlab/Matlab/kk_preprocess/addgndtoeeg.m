function out = addgndtoeeg(prefix, day, epochs, tets)
% Shantanu - Nov 2012.
% addgndtoeeg('bon', 1, 1:7, 1:30);

% DO ONE DAY AT A TIME - TO GET REFERENCE RIGHT
% (!!) Also, do for ALL TETRODES AT ONCE ..
        % sets eegstarttimes (and end sample) of all tetrodes to be the same (1 second + ceil of
        % first iterated tetrode, and -1 second a floor)
        
% Due to rounding, output eeggnd data length may be different in size from other tetrodes by a
% single sample.


if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epochs=1; %% Epochs
end
if nargin<4,
    tet=1; %
end


switch prefix
    case  'arn'        % Annabelle's Arnold
        directoryname = '/datatmp/kkay/Arn';
        reflist = 4  * ones(1,30);       % checked config file, all days use 4 as reference
    case 'bar'        % Annabelle's
        directoryname = '/data12/kkay/Bar';
        if day < 18
            reflist =  4 * ones(1,30);       % checked config file, days 1-17 use tet 4 as ref, while days 18-22 use tet 1
        else 
            reflist =  1 * ones(1,30);       % checked config file, days 1-17 use tet 4 as ref, while days 18-22 use tet 1
        end
    case 'cal'        % Annabelle's
        directoryname = '/data12/kkay/Cal';
        reflist = 16   * ones(1,30);       % checked config file, all days use 16 as reference    
    case 'dwi'        % Annabelle's
        directoryname = '/data12/kkay/Dwi';
        reflist = 7   * ones(1,30);       % checked config file, all days use 7 as reference    
        
        
    case 'Cor'
        directoryname = '/data12/kkay/Cor';
        reflist = 24 * ones(1,30);        % (!!) the indices in this vector should correspond to tetrode number        
    case 'cha'
        directoryname = '/data12/kkay/Cha/';
        reflist = 1 * ones(1,21);        % (!!) the indices in this vector should correspond to tetrode number
    case 'egy'
        directoryname = '/data12/mari/Egy/';
        reflist = 11 * ones(1,21);
    case 'bon'
        directoryname = '/data12/kkay/Bon';
        reflist = 30* ones(1,30);   % checked this for every day-epoch combo 5.2.13
    case 'fra'
        directoryname = '/vortexdata/kkay/Fra';
        reflist = 30 * ones(1,30);   % checked this for every day-epoch combo 5.2.13
    case 'HPa'
        rawdir = '/data25/sjadhav/HPExpt/HPa/';
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
        reflist = 3*ones(1,20);
    case 'HPb'
        rawdir = '/data25/sjadhav/HPExpt/HPb/';
        directoryname = '/data25/sjadhav/HPExpt/HPb_direct/';
        if (day==1) || (day==2)
            reflist = [7,7,7,7,7,7,7,11,11,11,11,11,11,11,7,7,7,7,7,7]; % dHp, PFC, iHp
        end
        if (day==3) || (day==4) || (day==5) || (day==6)
            reflist = [7,7,7,7,7,7,7,11,11,11,11,11,11,11,15,15,15,15,15,15];     % note here different references for different tetrodes
        end
        if (day==7) || (day==8)
            reflist = [7,7,7,7,7,7,7,13,13,13,13,11,13,13,15,15,15,15,15,15];
        end
             
end
dir2=directoryname;


%Adding Reference to EEG
%----------------------

cd([dir2,'/EEG/']);
refeeg=[]; curreeg=[];

flag_refchange=0;
% for d=1:length(days)
%     day=days(d);

if (day<10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end




for ep=1:length(epochs)
    startflag = 0;
    sampflag = 0;
    
    
            for t=1:length(tets)   % iterate through requested tetrodes
               
            tet=tets(t);            % select tetrode
            epoch=epochs(ep);       % select epoch
            eeggnd=[];          % initialize eeggnd (output variable)
            
            % retrieve ref tetrode # for this day
            ref = reflist(t);
            if (ref<10)
                refstring = ['0',num2str(ref)];
            else
                refstring = num2str(ref);
            end
            
            % take the reference tetrode to use for the new standard starttime and endtime for that epoch
            if startflag == 0
                reffilename = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',refstring];
                load(reffilename);
                reftimes = geteegtimes(eeg{day}{epoch}{ref});
                
                new_starttime=ceil(eeg{day}{epoch}{ref}.starttime + 0.5);   % designate new standard starttime
                new_endttime=floor(reftimes(end) - 0.5);                           % designate new standard enddtime
                
                ref_startsample=lookup(new_starttime,reftimes);
                ref_endsample=lookup(new_endttime,reftimes);
                refeeg = eeg{day}{epoch}{ref}.data(ref_startsample:ref_endsample);
                refsamprate=eeg{day}{epoch}{ref}.samprate;
                startflag = 1;
            end
           
            if ~ismember(tet,reflist)       % Only perform if Tet is NOT the reference / one of the references
            
            disp(['De-referencing Tet ',num2str(tet),', Epoch',num2str(epoch),'; Ref is ',num2str(ref)]);

            if (tet<10)
                currstring = ['0',num2str(tet)];
            else
                currstring = num2str(tet);
            end
            
            currname = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',currstring];
            load(currname);
            
            % retrieve data at new start and end times
            eegtimes=geteegtimes(eeg{day}{epoch}{tet});
            eeg_startsample=lookup(new_starttime,eegtimes);
            eeg_endsample=lookup(new_endttime,eegtimes);            
            trimmed_data = eeg{day}{epoch}{tet}.data(eeg_startsample:eeg_endsample);
            
            % quick fix correction code if lengths differ slightly
            if length(trimmed_data) ~= length(refeeg)
                if length(trimmed_data)-length(refeeg)==1
                   trimmed_data=trimmed_data(1:(end-1)); 
                elseif length(trimmed_data)-length(refeeg)==-1
                   refeeg=refeeg(1:(end-1));
                elseif length(trimmed_data)-length(refeeg)==2
                    trimmed_data=trimmed_data(1:(end-2));
                  disp(['Day ' num2str(day) ' 2 sample mismatch for Tet ',num2str(tet),' and Ref', num2str(ref)]);
                  disp(['length of trimmed data :' num2str(length(trimmed_data)) '  //  length of ref data :' num2str(length(refeeg))])                   
                elseif length(trimmed_data)-length(refeeg)==-2
                    refeeg=refeeg(1:(end-2));
                  disp(['Day ' num2str(day) ' 2 sample mismatch for Tet ',num2str(tet),' and Ref', num2str(ref)]);
                  disp(['length of trimmed data :' num2str(length(trimmed_data)) '  //  length of ref data :' num2str(length(refeeg))])
                else
                  disp(['Day ' num2str(day) ' >2 length mismatch for Tet ',num2str(tet),' and Ref', num2str(ref)]);
                  disp(['length of trimmed data :' num2str(length(trimmed_data)) '  //  length of ref data :' num2str(length(refeeg))])
                end
            end
            
            %in some animals, sampling rates with differ between dsps! in
            %this case use interpolation and let the tetrode w/ higher
            %samprate win
            samprate = eeg{day}{epoch}{tet}.samprate;
            diffsamprate = refsamprate-samprate;
            
            if (diffsamprate > 0) && abs(length(trimmed_data)-length(refeeg)) > 2
                trimmed_data = interp1(eegtimes(eeg_startsample:eeg_endsample),trimmed_data,...
                                       reftimes(ref_startsample:ref_endsample),'nearest');
                trimmed_data=trimmed_data';
            elseif (diffsamprate < 0) && abs(length(trimmed_data)-length(refeeg)) > 2
                if sampflag == 0
                    refeeg = interp1(reftimes(ref_startsample:ref_endsample),refeeg,...
                                       eegtimes(eeg_startsample:eeg_endsample),'nearest');
                    refeeg = refeeg';
                    sampflag = 1;
                end
            end
            
            % in some cases (older animals) you may encounter trimmed_data
            % and refeeg still having different #s of samples -- ignore in
            % this case
            diffsamples = length(trimmed_data)-length(refeeg);
            if diffsamples ~= 0
                disp(['Day ' num2str(day) ' -- unfixed --'  diffsamples ' sample mismatch for Tet ',num2str(tet),' and Ref', num2str(ref)]);
                disp(['    ... skipping ']);
                continue
            end
            
            
            eeggnd{day}{epoch}{tet}.data = trimmed_data + refeeg;
            eeggnd{day}{epoch}{tet}.starttime = new_starttime;
            eeggnd{day}{epoch}{tet}.samprate = eeg{day}{epoch}{tet}.samprate;
            eeggnd{day}{epoch}{tet}.depth = eeg{day}{epoch}{tet}.depth;
            eeggnd{day}{epoch}{tet}.fields = strvcat(eeg{day}{epoch}{tet}.fields,[' wrt GND']);
            eeggnd{day}{epoch}{tet}.descript = strvcat(eeg{day}{epoch}{tet}.descript,[' wrt GND']);
            % Now Save
            savename = [dir2,'/EEG/',prefix,'eeggnd',daystring,'-',num2str(epoch),'-',currstring];
            save(savename,'eeggnd');
            

            
            
        end % end epochs
        
    end % end NOT Ref
    
end % end tets

%end % end days




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



