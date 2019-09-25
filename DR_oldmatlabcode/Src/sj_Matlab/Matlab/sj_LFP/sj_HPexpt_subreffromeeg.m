function [ref] = sj_HPexpt_subreffromeeg(prefix, day, epochs, tets)
% Shantanu - Nov 2012.
% Subtract Ref EEG from Eeggnd, and save. From sj_HPexpt_addgndtoeeg.
% The eeggnd File Has to Exist for this to work

% sj_HPexpt_subreffromeeg('HPb', 1, 1:7, 8:14); % PFC tets wrt to Hip
% sj_HPexpt_subreffromeeg('HPa', 2, 1:5, 15:16); % PFC tets wrt to PFC

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
        ref = 19; % This is pfc reference. In this animal, all was recorded wrt to Hp Ref 3
        isrefwrtgnd = 'no'; % Was reference recorded wrt gnd originally? If not, then you need to load the gnd file 
        reftag = 'Ref19Pfc';
        reflist = 3; % Recorded wrt gnd
    case 'HPb'
        rawdir = '/data25/sjadhav/HPExpt/HPb/';
        directoryname = '/data25/sjadhav/HPExpt/HPb_direct/';
        ref = 7; % This is Hp Ref you want to subtract from
        isrefwrtgnd = 'yes'; % Was reference recorded wrt gnd originally? Then the original eeg file is already wrt gnd
        reftag = 'Ref7Hp';
        reflist = [7;11]; % Recorded wrt gnd
        
end
dir2=directoryname;


% Adding Reference to EEG
%----------------------

cd([dir2,'/EEG/']);
% for d=1:length(days)
%     day=days(d);

if (day<10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

for t=1:length(tets)
    tet=tets(t);
    if (ref<10)
        refstring = ['0',num2str(ref)];
    else
        refstring = num2str(ref);
    end

    for ep=1:length(epochs)
        
        epoch=epochs(ep);
        disp(['Doing Tet ',num2str(tet),', Epoch',num2str(epoch),'; Ref is ',num2str(ref)]);
        
        if strcmp(isrefwrtgnd,'yes') == 1
            % The eeg file already is wrt to gnd, so load it
            reffilename = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',refstring];
            load(reffilename);
            refeeg = eeg{day}{epoch}{ref}.data;
            clear eeg
            
        else
            % Load the eegGND file
            reffilename = [dir2,'/EEG/',prefix,'eeggnd',daystring,'-',num2str(epoch),'-',refstring];
            load(reffilename);
            refeeg = eeggnd{day}{epoch}{ref}.data;
            clear eeggnd
        end
            
        
        if (tet<10)
            currstring = ['0',num2str(tet)];
        else
            currstring = num2str(tet);
        end
        
          %     Check if Tet is not one of the original References
        if ~ismember(tet,reflist)
            currname = [dir2,'/EEG/',prefix,'eeggnd',daystring,'-',num2str(epoch),'-',currstring];
            load(currname);
            eeggnddata = eeggnd{day}{epoch}{tet}.data;
        else % If it was, eeggnd datais in the original file
            currname = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',currstring];
            load(currname);
            eeggnddata = eeg{day}{epoch}{tet}.data;
        end
        
        % Subtract Ref eeg from given eeg if lengths match
        if length(eeggnddata) == length(refeeg)
            eegref{day}{epoch}{tet}.data = eeggnddata - refeeg;
            % Now Save
            savename = [dir2,'/EEG/',prefix,'eeg',reftag,daystring,'-',num2str(epoch),'-',currstring]
            save(savename,'eegref');
        else
            disp(['Length mismatch for Tet ',num2str(tet),' and Ref', num2str(ref)]);
        end
        
        
    end % end epochs
    
    %     end % end NOT Ref
    
end % end tets

%end % end days




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



