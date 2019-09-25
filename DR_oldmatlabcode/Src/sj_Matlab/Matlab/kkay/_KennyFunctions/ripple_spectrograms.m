% given a selected tetrode (e.g. CA1), calculates spectrograms time locked
% to ripple (startind)

%% two ingredients: eeg structure (all epochs for a day, use Mattias'
%% loadeegstruct) + ripple structure
%% after rippledayprocess() and extractripples(),
%% obtain this struct: ripples{day}{epoch}{tetrode}.<field>

d=5;
windowsize_sec = 0.8;
Fs = 1500;
halfwindow_samp = (windowsize_sec*1500)/2;
ref_tetrode = 11;


% first collect all event times for day, from chosen reference tetrode

eventtimes=cell(1,6);
for e=1:6
    eventtimes{e}=[eventtimes{e}; ripples{d}{e}{ref_tetrode}.starttime]; % given a selected tetrode (e.g. CA1), calculates spectrograms time locked
end
   
% copy raw eeg into windows centered on event times

eegwindows=cell(1,7,14); % holds matrices of rows(events) and columns(times)

for e=1:6
    epoch_starttime=eeg{d}{e}{1}.starttime;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eeg{d}{e}{1}.data));
    for t=1:14
        for r=(1+5):(length(eventtimes{e})-5)          % iterate over all events except first and last 5
            centerindex=lookup(eventtimes{e}(r),epoch_times);
            eegwindows{1,e,t}=[eegwindows{1,e,t} ; ...
                                 eeg{d}{e}{t}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))'];

        end
    end
end

%%%%%%%%%%%% sharp-wave comparisons: plot raw eeg for example events across
%%%%%%%%%%%% various tetrodes


    e=6;  % set epoch manually
    timevec=((1:1201)-600)/1.5;
    for i=1:10
    figure
    hold on
    N=size(eegwindows{1,e,1},1);        % # events to choose from
    eventno=ceil(N*rand);
    string=sprintf('individual event no %d',eventno);
    title(string);
    for t=8:14
        if t==9          %% CA2
            plot(timevec,eegwindows{1,e,t}(eventno,:),'g','LineWidth',2);
        elseif t==8 || t==13             %% near CA2
            plot(timevec,eegwindows{1,e,t}(eventno,:),'LineWidth',2,'Color',[0 0.5 0]);
        elseif t==10             %% DG
             % plot(timevec,eegwindows{1,e,t}(eventno,:),'m','LineWidth',2);
        elseif t==12             %% CA3
             %plot(timevec,eegwindows{1,e,t}(eventno,:),'r','LineWidth',2);
        else                     %% CA1
            plot(timevec,eegwindows{1,e,t}(eventno,:),'k','LineWidth',2);
        string = sprintf('%d',t);
        title(string);
        end
    end
    end


%%%%%%%%%%%%%





% calculate spectrograms

allwindowedspectra=cell(1,7,14);          
movingwin = [100 10]/1000;            % chronux params
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];


for e=1:6
    for t=1:14
        flag=0;
        for r=(1+5):(length(eventtimes{e})-5)   %% ignore 5 events at the beginning and end
            [S,times,freqs,Serr]=mtspecgramc(eegwindows{1,e,t}(r-5,:),movingwin,params);
            
            if flag==0             % this if-clause initializes the 3D matrix of all spectrograms for a given tetrode
                allwindowedspectra{1,e,t}=nan(size(S,1),size(S,2),length(eventtimes{e})-10);
                flag=1;
            end
            
            allwindowedspectra{1,e,t}(:,:,r-5)=S;   % adds the S 2D matrix in the third dimension
        end
    end
end


% take mean event spectra (in epochs: meaneventspectra AND all epochs lumped: meaneventspectra_lump) 

meaneventspectra=cell(1,6,14);
for e=1:6
    for t=1:14
        meaneventspectra{1,e,t}=mean(allwindowedspectra{1,e,t},3);
    end
end

meaneventspectra_lump = cell(1,14);   % all epochs' events lumped% calculate spectrograms

allwindowedspectra=cell(1,7,14);          
movingwin = [100 10]/1000;            % chronux params
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];


for e=1:6
    for t=1:14
        flag=0;
        for r=(1+5):(length(eventtimes{e})-5)   %% ignore 5 events at the beginning and end
            [S,times,freqs,Serr]=mtspecgramc(eegwindows{1,e,t}(r-5,:),movingwin,params);
            
            if flag==0             % this if-clause initializes the 3D matrix of all spectrograms for a given tetrode
                allwindowedspectra{1,e,t}=nan(size(S,1),size(S,2),length(eventtimes{e})-10);
                flag=1;
            end
            
            allwindowedspectra{1,e,t}(:,:,r-5)=S;   % adds the S 2D matrix in the third dimension
        end
    end
end


% take mean event spectra (in epochs: meaneventspectra AND all epochs lumped: meaneventspectra_lump) 

meaneventspectra=cell(1,6,14);
for e=1:6
    for t=1:14
        meaneventspectra{1,e,t}=mean(allwindowedspectra{1,e,t},3);
    end
end

meaneventspectra_lump = cell(1,14);   % all epochs' events lumped
    for t=1:14
        dummy=[];
        for e=1:6
            dummy=cat(3,dummy,allwindowedspectra{1,e,t});   % concatenate events over epoch
        end
        meaneventspectra_lump{t}=mean(dummy,3);
    end

% find overall day mean and std spectra for each tetrode

meandayspectra=cell(1,14);
stddayspectra=cell(1,14);

for d=5
    for t=1:14
        dummy=[];
        for e=1:6
            [S_full,junkt,junkf,junkserr] = mtspecgramc(eeg{d}{e}{t}.data,movingwin,params);
            dummy=[dummy;S_full];        
        end
        meandayspectra{t}=mean(dummy,1);
        stddayspectra{t}=std(dummy,1);
    end
end

    for t=1:14
        dummy=[];
        for e=1:6
            dummy=cat(3,dummy,allwindowedspectra{1,e,t});   % concatenate events over epoch
        end
        meaneventspectra_lump{t}=mean(dummy,3);
    end


save('RIPPLE_SPECTROGRAMS_WORKSPACE_REF11','-v7.3');


% calculate ratio spectrogram for all lumped events from tetrode
%(using _lump)

ratiospectra=cell(1,14);
figure
for t=8:14
        ratiospectra{t}=bsxfun(@rdivide,meaneventspectra_lump{t},meandayspectra{t});
        subplot(1,14,t)
        imagesc(times,freqs,ratiospectra{t}',[0,4]);
            set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
end


% calculate z-scored spectra

zscorespectra_epochs=cell(1,6,14);
zscorespectra_lump=cell(1,14);

for e=6
    for t=1:14
        for r=1:size(allwindowedspectra{1,e,t},3)
        zscorespectra_epochs{1,e,t}(:,:,r)=bsxfun(@minus,allwindowedspectra{1,e,t}(:,:,r),meandayspectra{t});   
        zscorespectra_epochs{1,e,t}(:,:,r)=bsxfun(@rdivide,allwindowedspectra{1,e,t}(:,:,r),stddayspectra{t});
        end
    end
end

figure
title('z-scored spectrograms');
    for t=1:14              % lump
        dummy=[];
        for e=6
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,t});   % concatenate events over epoch
        end
        zscorespectra_lump{t}=mean(dummy,3);
                subplot(1,14,t)
        imagesc(times,freqs,zscorespectra_lump{t}',[0,5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
    end


    %%%%%%% manual plot of 3 tetrodes
        figure
            dummy=[];
        for e=1:6
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,8});   % concatenate events over epoch
        end
        zscorespectra_lump{1}=mean(dummy,3);
                subplot(1,3,1)
        imagesc(times,freqs,zscorespectra_lump{1}',[0,3.5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',7);
        title(string);
                dummy=[];
        for e=1:6
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,11});   % concatenate events over epoch
        end
        zscorespectra_lump{2}=mean(dummy,3);
                subplot(1,3,2)
        imagesc(times,freqs,zscorespectra_lump{2}',[0,3.5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',9);
        title(string);
                dummy=[];
        for e=1:6
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,14});   % concatenate events over epoch
        end
        zscorespectra_lump{t}=mean(dummy,3);
                subplot(1,3,3)
        imagesc(times,freqs,zscorespectra_lump{13}',[0,3.5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',13);
        title(string);
    %%%%%%%%%%%%%%%%
    
    
    
    
    % plot some random individual ripple-triggered z-scored spectrograms
    % from an epoch
    
    e=2;  % set epoch manually
    for i=1:1
    figure
    N=size(zscorespectra_epochs{1,e,t},3);        % # events to choose from
    eventno=ceil(N*rand);
    string=sprintf('individual event no %d',eventno);
    title(string);
    for t=1:14
        subplot(1,14,t)
        imagesc(times,freqs,zscorespectra_epochs{1,e,t}(:,:,eventno)',[0,5]);
        set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
        end
    end






