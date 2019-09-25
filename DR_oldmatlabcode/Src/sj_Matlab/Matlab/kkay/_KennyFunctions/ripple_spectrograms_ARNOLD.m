% given a selected tetrode (e.g. CA1), calculates spectrograms time locked
% to ripple (startind)

%% two ingredients: eeg structure (all epochs for a day, use Mattias'
%% loadeegstruct()   ) + ripple structure
%% after rippledayprocess() and extractripples(),
%% obtain this struct: ripples{day}{epoch}{tetrode}.<field>

windowsize_sec = 0.8;
Fs = 1500;
halfwindow_samp = (windowsize_sec*1500)/2;
ref_tetrode = 11;

% first collect all event times for day, from chosen reference tetrode

eventtimes=cell(1,4);
for e=1:4
    eventtimes{e}=[eventtimes{e}; ripples{5}{e}{ref_tetrode}.starttime];
end

% copy raw eeg into windows centered on event times

eegwindows=cell(5,4,14); % holds matrices of rows(events) and columns(times)

for e=1:4
    epoch_starttime=eeg{5}{e}{11}.starttime;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eeg{5}{e}{11}.data));
    for t=11
        for r=(1+3):(length(eventtimes{e})-3)          % iterate over events except last few
            centerindex=lookup(eventtimes{e}(r),epoch_times);
            eegwindows{5,e,t}=[eegwindows{5,e,t} ; ...
                                 eeg{5}{e}{t}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))'];

        end
    end
end


%%%%%%%%%%%% sharp-wave comparisons: plot raw eeg for example events across
%%%%%%%%%%%% various tetrodes


    e=1;  % set epoch manually
    for i=1:10
    figure
    hold on
    N=size(eegwindows{1,e,1},1);        % # events to choose from
    eventno=ceil(N*rand);
    string=sprintf('individual event no %d',eventno);
    title(string);
    for t=11
            plot(eegwindows{1,e,t}(eventno,:),'b');
        string = sprintf('%d',t);
        title(string);
        end
    end
    end


%%%%%%%%%%%%%





% calculate spectrograms

allwindowedspectra=cell(5,7,14);          
movingwin = [150 10]/1000;            % chronux params
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 350];
params.tapers = [3 5];


for e=1:4
    for t=11
        flag=0;
        for r=(1+3):(length(eventtimes{e})-3)   %% ignore 5 events at the beginning and end
            [S,times,freqs,Serr]=mtspecgramc(eegwindows{5,e,t}(r-3,:),movingwin,params);
            
            if flag==0             % this if-clause initializes the 3D matrix of all spectrograms for a given tetrode
                allwindowedspectra{5,e,t}=nan(size(S,1),size(S,2),length(eventtimes{e})-6);
                flag=1;
            end
            
            allwindowedspectra{5,e,t}(:,:,r-3)=S;   % adds the S 2D matrix in the third dimension
        end
    end
end


% take mean event spectra (in epochs: meaneventspectra AND all epochs lumped: meaneventspectra_lump) 

meaneventspectra=cell(1,7,14);
for e=1:4
    for t=11
        meaneventspectra{5,e,t}=mean(allwindowedspectra{5,e,t},3);
    end
end

meaneventspectra_lump = cell(1,14);   % all epochs' events lumped
    for t=11
        dummy=[]
        for e=1:4
            dummy=cat(3,dummy,allwindowedspectra{5,e,t});   % concatenate events over epoch
        end
        meaneventspectra_lump{t}=mean(dummy,3);
    end

% find overall day mean and std spectra for each tetrode

meandayspectra=cell(1,14);
stddayspectra=cell(1,14);

for d=5
    for t=11
        dummy=[]
        for e=1:4
            [S_full,junkt,junkf,junkserr] = mtspecgramc(eeg{5}{e}{t}.data,movingwin,params);
            dummy=[dummy;S_full];        
        end
        meandayspectra{t}=mean(dummy,1);
        stddayspectra{t}=std(dummy,1);
    end
end





% calculate z-scored spectra

zscorespectra_epochs=cell(5,7,14);
zscorespectra_lump=cell(5,14);

for e=1:4
    for t=11
        for r=1:size(allwindowedspectra{5,e,t},3)
        zscorespectra_epochs{5,e,t}(:,:,r)=bsxfun(@minus,allwindowedspectra{5,e,t}(:,:,r),meandayspectra{t});   
        zscorespectra_epochs{5,e,t}(:,:,r)=bsxfun(@rdivide,allwindowedspectra{5,e,t}(:,:,r),stddayspectra{t});
        end
    end
end

figure
title('z-scored spectrograms');
    for t=11              % lump
        dummy=[];
        for e=1:4
            dummy=cat(3,dummy,zscorespectra_epochs{5,e,t});   % concatenate events over epoch
        end
        zscorespectra_lump{t}=mean(dummy,3);
        imagesc(times,freqs,zscorespectra_lump{t}',[0,1]);
            set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
    end

    
    