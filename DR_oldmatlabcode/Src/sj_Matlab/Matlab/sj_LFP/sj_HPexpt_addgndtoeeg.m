function [ref] = sj_HPexpt_addgndtoeeg(prefix, day, epochs, tets)
% DO SJ_EEGUPDATE FOR TETS ON NEW DSPS (USUALLY 15 ONWARDS) TO MATCH LENGHTS BEFORE RUNNING
% Shantanu - Nov 2012.
% sj_HPexpt_addgndtoeeg('HPa', 8, 1:5, 17:20);
% sj_HPexpt_addgndtoeeg('HPb', 1, 1:7, 1:20);
% sj_HPexpt_addgndtoeeg('HPb', 2, 1:5, 15:20);

% DO ONE DAY AT A TIME - TO GET REFERENCE RIGHT

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
            reflist = [7,7,7,7,7,7,7,11,11,11,11,11,11,11,15,15,15,15,15,15];
        end
        if (day==7) || (day==8)
            reflist = [7,7,7,7,7,7,7,13,13,13,13,11,13,13,15,15,15,15,15,15];
        end
    case 'HPc'
        rawdir = '/data25/sjadhav/HPExpt/HPc/';
        directoryname = '/data25/sjadhav/HPExpt/HPc_direct/';
        if (day==1) || (day==2) || (day==3) || (day==4)
            reflist = [7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,17,17,17,17,17,17]; % dHp, PFC, iHp
        end
        if (day==5) || (day==6) || (day==7) || (day==8)
            reflist = [7,7,7,7,7,7,7,13,13,13,13,13,13,13,13,13,17,17,17,17,17,17]; % dHp, PFC, iHp
        end
end
dir2=directoryname;

disp(['Adding Ground ... ']);

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

for t=1:length(tets)
    tet=tets(t);
    ref = reflist(tet);   % V. IMP - SHOULD BE "TET" IN BRACKET
    if (ref<10)
        refstring = ['0',num2str(ref)];
    else
        refstring = num2str(ref);
    end
    
    % Only do if Tet is not one of the References, else the eeg is already the Gnd
    % --------------------------------------------
    if ~ismember(tet,reflist)
        
        %             if t==1
        %                 ref = reflist(t); % For first tet
        %                 flag_refchange=1;
        %             end
        %             % For 2nd tet onwards, see whether this is a new reference from previous.
        %             if t > 1
        %                 currref = reflist(t);
        %                 if currref ~= ref
        %                     ref = currref; % update and load new if it has changed
        %                     flag_refchange=1;
        %                 else
        %                     flag_refchange=0;
        %                 end
        %             end
        %             if flag_refchange ==1
        %                 if (ref<10)
        %                     refstring = ['0',num2str(ref)];
        %                 else
        %                     refstring = num2str(ref);
        %                 end
        %                 reffilename = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',refstring];
        %                 load(reffilename);
        %                 refeeg = eeg{day}{epoch}{ref}.data;
        %             end
        
        
        for ep=1:length(epochs)
            
            % Load and initialize eeggnd
            % -------------------------
            eeggnd=[]; 
            epoch=epochs(ep);
            disp(['Doing Tet ',num2str(tet),', Epoch',num2str(epoch),'; Ref is ',num2str(ref)]);
            reffilename = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',refstring];
            load(reffilename);
            refeeg = eeg{day}{epoch}{ref}.data;
            
            % Load curr tet
            % ----------------
            if (tet<10)
                currstring = ['0',num2str(tet)];
            else
                currstring = num2str(tet);
            end
            currname = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',currstring];
            load(currname);
            
            % Add back Ref data if length match
            % -----------------------------------
            
            if length(eeg{day}{epoch}{tet}.data) == length(refeeg)
                eeggnd{day}{epoch}{tet}.data = eeg{day}{epoch}{tet}.data + refeeg;
                eeggnd{day}{epoch}{tet}.starttime = eeg{day}{epoch}{tet}.starttime;
                eeggnd{day}{epoch}{tet}.samprate = eeg{day}{epoch}{tet}.samprate;
                eeggnd{day}{epoch}{tet}.depth = eeg{day}{epoch}{tet}.depth;
                eeggnd{day}{epoch}{tet}.fields = strvcat(eeg{day}{epoch}{tet}.fields,[' wrt GND']);
                eeggnd{day}{epoch}{tet}.descript = strvcat(eeg{day}{epoch}{tet}.descript,[' wrt GND']);
                % Now Save
                % ----------
                savename = [dir2,'/EEG/',prefix,'eeggnd',daystring,'-',num2str(epoch),'-',currstring];
                save(savename,'eeggnd');
            else
                disp(['Length mismatch for Tet ',num2str(tet),' and Ref', num2str(ref)]);
            end
            
            
        end % end epochs
        
    else %(if tet is a reference)
        
        disp(['Tet ',num2str(tet),' is a reference. So copying the same eeg to the Gnd file']);
        
        for ep=1:length(epochs)
            eeggnd=[];
            epoch=epochs(ep);
            disp(['Copying Ref Tet ',num2str(tet),', Epoch',num2str(epoch)]);
            if (tet<10)
                currstring = ['0',num2str(tet)];
            else
                currstring = num2str(tet);
            end
            currname = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',currstring];
            load(currname);
            eeggnd{day}{epoch}{tet} = eeg{day}{epoch}{tet};
            savename = [dir2,'/EEG/',prefix,'eeggnd',daystring,'-',num2str(epoch),'-',currstring];
            save(savename,'eeggnd');
        end % end epochs

        
    end % end NOT Ref
    
end % end tets

%end % end days




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



