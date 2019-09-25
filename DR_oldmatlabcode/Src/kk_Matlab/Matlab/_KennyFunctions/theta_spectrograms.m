% given a selected tetrode (e.g. CA1), calculates spectrograms time locked
% to hightheta (startind)

%% two ingredients: eeg structure (all epochs for a day, use Mattias'
%% consolidator, e.g. "eeg_day_1.mat") + theta structure
%% after thetadayprocess() and extracthightheta(),
%% obtain this struct: hightheta{day}{epoch}{tetrode}.<field>

windowsize_sec = 0.8;
Fs = 1500;
halfwindow_samp = (windowsize_sec*1500)/2;
ref_tetrode = 7;

% first collect all event times for day, from chosen reference tetrode

eventtimes=cell(1,7);
for e=1:7
    eventtimes{e}=[eventtimes{e}; hightheta{1}{e}{ref_tetrode}.starttime]
end

% copy raw eeg into windows centered on event times

eegwindows=cell(1,7,14); % holds matrices of rows(events) and columns(times)

for e=1:7
    epoch_starttime=eeg{1}{e}{1}.starttime;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eeg{1}{e}{1}.data));
    for t=1:14
        for r=(1+1):(length(eventtimes{e})-1)          % iterate over all events
            centerindex=lookup(eventtimes{e}(r),epoch_times);
            eegwindows{1,e,t}=[eegwindows{1,e,t} ; ...
                                 eeg{1}{e}{t}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))'];

        end
    end
end

%%%%%%%%%%%% raw LFP comparisons: plot raw eeg for example events across
%%%%%%%%%%%% various tetrodes

    e=1;  % set epoch manually
    timevec=((1:1201)-600)/1.5;
    N=size(eegwindows{1,e,1},1);        % # events to choose from
    for r=1:N
        figure
        hold on
    string=sprintf('individual event no %d',r);
    title(string);
    for t=[1 2 7 9 13]
        if t==13
            plot(timevec,eegwindows{1,e,t}(r,:),'r');
        else
            plot(timevec,eegwindows{1,e,t}(r,:),'b');
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


for e=1:7
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


% take mean event spectra (in epochs: meaneventspectra OR all epochs lumped: meaneventspectra_lump) 

meaneventspectra=cell(1,7,14);
for e=1:7
    for t=1:14
        meaneventspectra{1,e,t}=mean(allwindowedspectra{1,e,t},3);
    end
end

meaneventspectra_lump = cell(1,14);   % all epochs' events lumped
    for t=1:14
        dummy=[]
        for e=1:7
            dummy=cat(3,dummy,allwindowedspectra{1,e,t});   % concatenate events over epoch
        end
        meaneventspectra_lump{t}=mean(dummy,3);
    end

% find overall day mean and std spectra for each tetrode

meandayspectra=cell(1,14);
stddayspectra=cell(1,14);

for d=1
    for t=1:14
        dummy=[]
        for e=1:7
            [S_full,junkt,junkf,junkserr] = mtspecgramc(eeg{1}{e}{t}.data,movingwin,params);
            dummy=[dummy;S_full];        
        end
        meandayspectra{t}=mean(dummy,1);
        stddayspectra{t}=std(dummy,1);
    end
end



% calculate z-scored spectra

zscorespectra_epochs=cell(1,7,14);
zscorespectra_lump=cell(1,14);

for e=[1 2 6 7]           %% ** no events on day 1 for epochs 3 through 5!
    for t=1:14
        for r=(1+1):(size(allwindowedspectra{1,e,t},3)-1)
        zscorespectra_epochs{1,e,t}(:,:,r)=bsxfun(@minus,allwindowedspectra{1,e,t}(:,:,r),meandayspectra{t});   
        zscorespectra_epochs{1,e,t}(:,:,r)=bsxfun(@rdivide,allwindowedspectra{1,e,t}(:,:,r),stddayspectra{t});
        end
    end
end

figure
title('z-scored spectrograms');
    for t=1:14              % lump
        dummy=[];
        for e=1:7
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,t});   % concatenate events over epoch
        end
        zscorespectra_lump{t}=mean(dummy,3);
                subplot(1,14,t)
        imagesc(times,freqs,zscorespectra_lump{t}',[0,1]);
            set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
    end
    
    
    
    %%%%%%%%%%% for theta, narrow focus on 0-15 Hz band
    figure
title('z-scored spectrograms');

        dummy2=[];
for t=1:14              % lump
        dummy=[];
        for e=1:7
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,t});   % concatenate events over epoch
        end
        zscorespectra_lump{t}=mean(dummy,3);
        dummy2=[dummy2;zscorespectra_lump{t}];
                subplot(14,1,t)
        imagesc(times,freqs(1:5),zscorespectra_lump{t}(:,1:5)',[0,1]);
            set(gca,'YDir','normal');
        title(string);
 end
    figure
        imagesc(dummy2',[0,1]);
    

    
    

    %%%%%%% manual plot of 3 tetrodes
        figure
            dummy=[];
        for e=1:7
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,1});   % concatenate events over epoch
        end
        zscorespectra_lump{1}=mean(dummy,3);
                subplot(1,3,1)
        imagesc(times,freqs,zscorespectra_lump{1}',[0,1.5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',1);
        title(string);
                dummy=[];
        for e=1:7
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,7});   % concatenate events over epoch
        end
        zscorespectra_lump{2}=mean(dummy,3);
                subplot(1,3,2)
        imagesc(times,freqs,zscorespectra_lump{2}',[0,1.5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',3);
        title(string);
                dummy=[];
        for e=1:7
            dummy=cat(3,dummy,zscorespectra_epochs{1,e,13});   % concatenate events over epoch
        end
        zscorespectra_lump{t}=mean(dummy,3);
                subplot(1,3,3)
        imagesc(times,freqs,zscorespectra_lump{13}',[0,1.5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',13);
        title(string);
    %%%%%%%%%%%%%%%%
    
    
    
    
    % plot some random individual theta-triggered z-scored spectrograms
    % from an epoch
    
    e=7;  % set epoch manually
    for i=1:10
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






