function [params] = sj_HPexpt_baselinespecgram(prefix, days, epochs, tets, do_wrtgnd,varargin)
% Shantanu - Nov 2012. 
% From Kenny and Maggies - event_spectrograms .m and calcriptriggerredspectrograms.m respectively
% This is for getting baseline values for given eeg tets for normalization when needed. Save These. 

% sj_HPexpt_baselinespecgram('HPa', 8, [1:5], [1:20], 0);
% sj_HPexpt_baselinespecgram('HPa', 2, [1:5], [1:20], 0);

% sj_HPexpt_baselinespecgram('HPb', 1, [1:7], 1:20, 0)
% sj_HPexpt_baselinespecgram('HPa', 2, [1:5], [14,15], 0);
% sj_HPexpt_baselinespecgram('HPb', 2, [1:5], [4,9,16], 0);
% sj_HPexpt_baselinespecgram('HPb', 6, [1:5], [4,9,16], 0);



if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epoch=1; %% Epochs 
end
if nargin<4,
    tets=1; %
end
if nargin<5,
    do_wrtgnd=0; % Whether to also do with respect to ground
end


% Define Chronux params
% -------------------------------------------
movingwin = [1000 100]/1000; %movingwin = [100 10]/1000;                
params.Fs = 1500;
params.fpass = [0 40]; % params.fpass = [0 400];
params.tapers = [3 5];
params.err = [2 0.05];

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'movingwin'
            movingwin = varargin{option+1};
        case 'fpass'
    	    params.fpass = varargin{option+1};
    end
end

savetag = 'tmp';

if params.fpass(2) == 400
    savetag = '';
    %movingwin = [100 10]/1000; 
end
if params.fpass(2) == 100
    savetag = 'mid';
    %movingwin = [400 40]/1000;
end
if params.fpass(2) == 40
    %movingwin = [1000 100]/1000;
    savetag = 'low';
end
if params.fpass(2) <= 10
    %movingwin = [8000 800]/1000; 
    savetag = 'floor';
   
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
end
dir2=directoryname;

% % Load times file - if it will be needed
% % ---------------------------------------
% currdir = pwd;
% cd(rawdir);
% dayfolders = dir;
% daystr = sprintf('%02d', day);
% for i = 3:length(dayfolders)
%     if dayfolders(i).isdir
%         if strcmp(daystr,dayfolders(i).name(1:2))
%             disp(upper(dayfolders(i).name))
%             cd(dayfolders(i).name);
%             load times;
%         end
%     end
% end
% cd(currdir);
% Now Getting Range directly from times file
% userange = ranges(epoch+1,:); % 1st row is allepochs. So you need +1
% usename = names{epoch+1}(end-15:end);

% Get EEGs and spectrogram it
% ------------------------------
cd([directoryname,'/EEG/']);
savedir = [directoryname,'/EEGSpec/'];

for d=1:length(days)
    
    day = days(d);
    if (day<10)
        daystring = ['0',num2str(day)];
    else
        daystring = num2str(day);
    end
    
    
    for t=1:length(tets)
        tet=tets(t);
        
        if (tet<10)
            tetstring = ['0',num2str(tet)];
        else
            tetstring = num2str(tet);
        end
        eegspec = []; dummy=[]; % Save 1 file for each tet in a day 
        if do_wrtgnd==1
            eeggndspec = []; dummy_gnd = [];
        end
         
        flaggnd=0;
        for ep=1:length(epochs)
            eeg=[]; eeggnd=[];
            epoch=epochs(ep);
            disp(['Doing Day ',num2str(day) ', Ep ',num2str(epoch),', Tet ',num2str(tet)]);
            curreegfile = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',tetstring];
            load(curreegfile);
            %eegspec{day}{epoch}{tet}.specgram = mtspecgramc(eeg{day}{epoch}{tet}.data,movingwin,params);
            specgram = mtspecgramc(eeg{day}{epoch}{tet}.data,movingwin,params);
            % Dont save specgram in file - its too big
            % --------------------------------------------
            eegspec{day}{epoch}{tet}.meanspec = mean(specgram,1);
            eegspec{day}{epoch}{tet}.stdspec = std(specgram,1);
            dummy=[dummy; specgram]; % For combining across epochs            
                       
            % SHIFT GND TO ANOTHER FILE TO AVOID CONFUSION
%             % If do the gnd channel also          
%             if do_wrtgnd==1
%                 curreeggndfile = [dir2,'/EEG/',prefix,'eeggnd',daystring,'-',num2str(epoch),'-',tetstring];
%                 if (exist([curreeggndfile,'.mat'],'file'))==2
%                     flaggnd=1;
%                     load(curreeggndfile);
%                     specgram_gnd = mtspecgramc(eeggnd{day}{epoch}{tet}.data,movingwin,params);
%                     eeggndspec{day}{epoch}{tet}.meanspec = mean(specgram_gnd,1);
%                     eeggndspec{day}{epoch}{tet}.stdspec = std(specgram_gnd,1);
%                     dummy_gnd=[dummy_gnd;specgram_gnd];
%                     % Dont save specgram in file - its too big
%                     %eeggndspec{day}{epoch}{tet}.specgram  = [];
%                 end
%             end
            
            % Dont save specgram in file - its too big
            % --------------------------------------------
            %eegspec{day}{epoch}{tet}.specgram  = [];
            
        end % end epochs
        % Also save mean and std for whole day. Save in fields for 1st epoch - [FIXED This has become last epoch now - by accident. Fix later]
        eegspec{day}{1}{tet}.meandayspec=mean(dummy,1);
        eegspec{day}{1}{tet}.stddayspec=std(dummy,1);
        % Save File for current day and tet
        savefile = [savedir,prefix,'eegspec',savetag,daystring,'-Tet',tetstring]
        save(savefile,'eegspec');
%         if ((do_wrtgnd==1) && (flaggnd==1))
%             eeggndspec{day}{1}{tet}.meandayspec=mean(dummy_gnd,1);
%             eeggndspec{day}{1}{tet}.stddayspec=std(dummy_gnd,1);
%             %eegspec=eeggndspec; % To have the same name - Do this
%             savefile = [savedir,prefix,'eeggndspec',daystring,'-Tet',tetstring]
%             save(savefile,'eeggndspec');
%         end
        
    end % end tets
    

    
end % end days


















i=1;
%set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



