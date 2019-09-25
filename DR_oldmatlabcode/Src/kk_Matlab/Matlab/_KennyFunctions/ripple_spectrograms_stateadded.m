% given a selected tetrode (e.g. CA1), calculates spectrograms time locked
% to ripple (startind)

%% THREE ingredients:
  %%  1. eeg struct     (all epochs for a day, use Mattias' loadeegstruct)
  %%  2. ripple struct  (after rippledayprocess and rippleextract)
  %%  3. pos struct     (NOT rawpos)

d=4;                         %%% currently set for Bashir day 04 preliminary analysis
windowsize_sec = 0.8;
Fs = 1500;
halfwindow_samp = (windowsize_sec*1500)/2;
ref_tetrode = 9;                 


% first collect all event times for day, from chosen reference tetrode

eventtimes=cell(1,7);
for e=2:3
    eventtimes{e}=[eventtimes{e}; ripples{d}{e}{ref_tetrode}.starttime]     % given a selected tetrode (e.g. CA1), calculates spectrograms time locked
end
   
% copy raw eeg into windows centered on event times

eegwindows=cell(1,7,14,3); % holds matrices of rows(events) and columns(times)  *** state dimension added

for e=2:3
    epoch_starttime=eeg{d}{e}{1}.starttime;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eeg{d}{e}{1}.data));
    for t=8:14
        for r=(1+5):(length(eventtimes{e})-5)          % iterate over all events except first and last 5
            centerindex=lookup(eventtimes{e}(r),epoch_times);
            velocityindex=lookup(eventtimes{e}(r),pos{d}{e}.data(:,1));            %%%%%%%%%%%%%%%   added for state
            velocity=pos{d}{e}.data(velocityindex,5);
            s=0;
            if velocity<2           %% state 1    %% immobile state
                s=1;
            elseif velocity>8       %% state 2   %% run state
                s=2;
            else                    %% state 3   %% intermediate state
                s=3;
            end
            eegwindows{1,e,t,s}=cat(1,eegwindows{1,e,t,s}, ...
                eeg{d}{e}{t}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))');
        end
    end
end

%%% count states


statecounts=zeros(3,7);
for e=2:3
        for r=(1+5):(length(eventtimes{e})-5)
            velocityindex=lookup(eventtimes{e}(r),pos{d}{e}.data(:,1));            %%%%%%%%%%%%%%%   added for state
            velocity=pos{d}{e}.data(velocityindex,5);
            if velocity<2           %% state 1
                statecounts(1,e)=statecounts(1,e)+1;
            elseif velocity>8       %% state 2
                statecounts(2,e)=statecounts(2,e)+1;
            else                    %% state 3
                statecounts(3,e)=statecounts(3,e)+1;
            end
        end
end

            
 %%%%%%%%%%%%%%%%

% 
%%%%%%%%%%%% sharp-wave comparisons: plot raw eeg for example events across
%%%%%%%%%%%% various tetrodes


    e=3;  % set epoch manually
    timevec=((1:1201)-600)/1.5;
    for i=1:5
    figure
    hold on
    N=size(eegwindows{1,e,1,1},1);        % # events to choose from
    eventno=ceil(N*rand);
    string=sprintf('individual event no %d',eventno);
    title(string);
    for t=[8 9 13 11 14 12]
        if t==8 || t==9 || t==13
            plot(timevec,eegwindows{1,e,t,1}(eventno,:),'r');
        else
            plot(timevec,eegwindows{1,e,t,1}(eventno,:),'b');
        string = sprintf('%d',t);
        title(string);
        end
    end
    end

% 
% %%%%%%%%%%%%%





% calculate spectrograms

allwindowedspectra=cell(1,7,14,3);          
movingwin = [100 10]/1000;            % chronux params
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];


for e=2:3
    for t=8:14
        for s=1:3
            flag=0;
            for r=1:size(eegwindows{1,e,t,s},1);
                [S,times,freqs,Serr]=mtspecgramc(eegwindows{1,e,t,s}(r,:),movingwin,params);
                if flag==0             % this if-clause initializes the 3D matrix of all spectrograms for a given tetrode
                    allwindowedspectra{1,e,t,s}=nan(size(S,1),size(S,2),size(eegwindows{1,e,t,s},1));
                    flag=1;
                end
                allwindowedspectra{1,e,t,s}(:,:,r)=S;   % adds the S 2D matrix in the third dimension
            end
        end
    end
end

disp('hello')

% take mean event spectra, while segregating state

meaneventspectra=cell(1,7,14);
for e=2:3
    for t=8:14
        for s=1:3
            meaneventspectra{1,e,t,s}=mean(allwindowedspectra{1,e,t,s},3);       
        end
    end
end

meaneventspectra_statelump = cell(1,14,3);   % OVERALL mean event spectra
    for t=8:14
        for s=1:3
             dummy=[];
             for e=2:3
                 dummy=cat(3,dummy,allwindowedspectra{1,e,t,s});   % concatenate events over epoch
             end
             meaneventspectra_statelump{1,t,s}=mean(dummy,3);
        end
   end

% find overall day mean and std spectra, over all windows, for each tetrode

meandayspectra=cell(1,14);
stddayspectra=cell(1,14);

for d=4
    for t=8:14
        dummy=[]
        for e=2:3
            [S_full,junkt,junkf,junkserr] = mtspecgramc(eeg{4}{e}{t}.data,movingwin,params);
            dummy=[dummy;S_full];        
        end
        meandayspectra{t}=mean(dummy,1);
        stddayspectra{t}=std(dummy,1);
    end
end

save('Bashir_4.25.12_state_rippletriggeredspectrograms_ref9','-v7.3');


% % 
% % % calculate ratio spectrogram for all lumped events from tetrode
% % %(using _lump)
% % 
% ratiospectra=cell(1,14);
% figure
% for t=[8 9 13 11 14 12]
%         ratiospectra{t}=bsxfun(@rdivide,meaneventspectra_statelump{1,t,2},meandayspectra{t});
%         subplot(1,14,t)
%         imagesc(times,freqs,ratiospectra{t}',[-0.2,5]);
%             set(gca,'YDir','normal');
%         string = sprintf('%d',t);
%         title(string);
% end


% calculate z-scored spectra

zscorespectra_epochs=cell(1,7,14,3);

for e=2:3
    for t=8:14
        for s=1:3
            for r=1:size(allwindowedspectra{1,e,t,s},3)
                zscorespectra_epochs{1,e,t,s}(:,:,r)=bsxfun(@minus,allwindowedspectra{1,e,t,s}(:,:,r),meandayspectra{t});   
                zscorespectra_epochs{1,e,t,s}(:,:,r)=bsxfun(@rdivide,allwindowedspectra{1,e,t,s}(:,:,r),stddayspectra{t});
            end
        end
    end
end

%% zscorespectra_lumpmean

zscorespectra_lumpmean=cell(1,14,3);     %% no epochs, maintains state
zscorespectra_lumpmean2=cell(1,14);      %% discards state

    for t=8:14              % lump
        dummy2=[];
        for s=1:3
        dummy=[];
        for e=2:3
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,t,s});   % concatenate events over epoch
        end
        zscorespectra_lumpmean{1,t,s}=mean(dummy,3);
        zscorespectra_lumpmean2{1,t}=cat(3,dummy,zscorespectra_lumpmean2{1,t});
        end
    end
    
    
% %%%%%%%%%% plot z-score spectra, all
%     figure
% for t=8:14
%       subplot(1,14,t)
%         imagesc(times,freqs,mean(zscorespectra_lumpmean2{t},3)',[-0.2,5]);
%             set(gca,'YDir','normal');
%         string = sprintf('%d',t);
%         title(string);
% end
% colorbar


%% plot z-score spectra, each state
for s=1:3
figure
string = sprintf('state %d',s);
title(string);
for t=8:14
        subplot(1,14,t)
        imagesc(times,freqs,zscorespectra_lumpmean{1,t,s}',[-0.2,5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
end
colorbar
end







        
        
        
        
        
        
        
        
        
        
        
        
    
    % plot some random individual ripple-triggered z-scored spectrograms
    % from an epoch
    
    e=2;  % set epoch manually
    for i=1:5
    figure
    N=size(zscorespectra_epochs{1,e,t},3);        % # events to choose from
    eventno=ceil(N*rand);
    string=sprintf('individual event no %d',eventno);
    title(string);
    for t=8:14
        subplot(1,14,t)
        imagesc(times,freqs,zscorespectra_epochs{1,e,t}(:,:,eventno)');
        set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
        end
    end






