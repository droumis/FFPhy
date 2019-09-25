global subplot_count;
subplot_count = 1;

Veqn = '>=0';
minV =  str2num(Veqn(end));
maxstage = 3; % [1 2 3]
minVPF = 0; %cm/sec
minPeakPF = 3;
lessthan=0;
includestates = 6;

%Animal selection
%-----------------------------------------------------
animals = {'N2'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filterionno

days='[1]';%,'1:10';
%days = '[1:1]';
%days = '[9:9]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = [''];

%epochtype='Run';

%epochfilter{1} = ['isequal($epochtype, ''Run'')'];
epochfilter{1} = ['isequal($epoch,  4)'];

%cellfilter = '(isequal($area, ''ACC'') && ($meanrate >0 ) )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0 ) )'  ;

%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

%timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%    {'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))},{'JY_getbarrier','($barrier== 0)'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { 'JY_getlinvelocity', '$velocity > 0' };
tetrodefilter=['isequal($area, ''HP'')'];
%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);

f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'eegtetrodes',tetrodefilter);

%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'JY_singleepochanal';

f = setfilteriterator(f,iterator);

out=setfilterfunction(f, 'JY_trialeeg', {'data','meaneegspectrograms','ripple'},'ripplereference',16);

out=runfilter(out);

outdata=out.output{1,1};

% process each epoch separately

% based on event_spectrograms.m from KKay

%% plot speed and power spectrum

% % generate index to show day epoch and tetrode channel
% trialno=9;
% 
% %figure;
% tetrodes=[15:21];
% 
% speed=out.output{1,1}{1,1}.trialspeed{trialno};
% 
% figure;
% set(gcf,'position',[0 0 600 600]);
% set(gcf,'PaperPositionMode','auto');
% 
% nrows=size(tetrodes,2)+1;
% columns=1;
% %gap_h=0.0000000001;
% %gap_w=0.00000001;
% gap_h=0.005;
% gap_w=0.05;
% 
% marg_h=[0.05 0.10];
% marg_w=[0.1 0.1];
% ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
% 
% axes(ha(1));
% 
% plot(speed);
% xlim([0 length(speed)]);
% 
% ax=1;
% axes(ha(2));
% 
% plot(speed);
% xlim([0 length(speed)]);
% %singleeeg=out.output{1,1}{1,1}.intertrialeeg{1,trialno}(channel,:);
% %speed=out.output{1,1}{1,1}.intertrialspeed{trialno};
% 
% for ii=2:nrows;
%     for channel=tetrodes(ii-1);
%         ax=ii;
%         axes(ha(ax));
%         
%         singleeeg=out.output{1,1}{1,1}.eeg{1,trialno}(channel,:);
%         refeeg=out.output{1,1}{1,1}.raweeg{channel};
%         
%         movingwin = [200 20]/1000;
%         params.Fs = 1500;
%         params.err = [2 0.05];
%         params.fpass = [75 350];
%         params.tapers = [3 5];
%         
%         % calculate baseline epoch eeg
%         
%         [S,t,freqs]=mtspecgramc(refeeg,movingwin,params);
%         Smean=mean(S); Sstd=std(S);
%         [S2,t,freqs] = mtspecgramc(singleeeg,movingwin,params);
%         S2=bsxfun(@minus,S2,Smean);
%         S2=bsxfun(@rdivide,S2,Sstd);
%         
%         imagesc(t,freqs,S2',[0 0.5]);
%         %imagesc(tmt-win(1),fmt,epocheegmean',[0, 20]);
%         set(gca,'YDir','normal');
%         %string = sprintf('%d',tmt);
%         ylabel(sprintf('Tetrode %s trial %s',num2str(channel),num2str(trialno)));
%         colormap('jet')
%         %colorbar
%         set(gca,'xtick',[]);
%     end
% end


%% plot raw eeg trace intertrial

% generate index to show day epoch and tetrode channel
for trialno=6:9
%figure;
tetrodes=[1:21];
speed=out.output{1,1}{1,1}.intertrialspeed{trialno};
ripple=out.output{1,1}{1,1}.intertrialripple{trialno};

figure; set(gcf,'position',[0 0 600 600]); set(gcf,'PaperPositionMode','auto');

nrows=size(tetrodes,2)+1;
columns=1;
%gap_h=0.0000000001;
%gap_w=0.00000001;
gap_h=0.005; gap_w=0.05;
marg_h=[0.05 0.10]; marg_w=[0.1 0.1];

ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
axes(ha(1));
plot(speed);
xlim([0 length(speed)]); title(sprintf('Trial %s',num2str(trialno)));

axes(ha(2));
plot(ripple)
xlim([0 length(ripple)]);
%singleeeg=out.output{1,1}{1,1}.intertrialeeg{1,trialno}(channel,:);
%speed=out.output{1,1}{1,1}.intertrialspeed{trialno};

for ii=3:nrows;
    for channel=tetrodes(ii-1);
        ax=ii;
        axes(ha(ax));
        
        singleeeg=out.output{1,1}{1,1}.intertrialeeg{1,trialno}(channel,:);

        plot(singleeeg)
        %ylim([min(singleeeg) max(singleeeg)]);
         ylim([-400 400]);
         xlim([0 length(singleeeg)]);
        ylabel(sprintf('Tetrode %s',num2str(channel),num2str(trialno)));
        set(gca,'xtick',[],'ytick',[]);
    end
end
end

%% plot raw eeg trace trial

% generate index to show day epoch and tetrode channel
for trialno=1:6
%figure;
tetrodes=[1:21];
speed=out.output{1,1}{1,1}.trialspeed{trialno};
ripple=out.output{1,1}{1,1}.trialripple{trialno};

figure; set(gcf,'position',[0 0 600 600]); set(gcf,'PaperPositionMode','auto');

nrows=size(tetrodes,2)+1;
columns=1;
%gap_h=0.0000000001;
%gap_w=0.00000001;
gap_h=0.005; gap_w=0.05;
marg_h=[0.05 0.10]; marg_w=[0.1 0.1];

ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
axes(ha(1));
plot(speed);
xlim([0 length(speed)]); title(sprintf('Trial %s',num2str(trialno)));

axes(ha(2));
plot(ripple)
xlim([0 length(ripple)]);
%singleeeg=out.output{1,1}{1,1}.intertrialeeg{1,trialno}(channel,:);
%speed=out.output{1,1}{1,1}.intertrialspeed{trialno};

for ii=3:nrows;
    for channel=tetrodes(ii-1);
        ax=ii;
        axes(ha(ax));
        
        singleeeg=out.output{1,1}{1,1}.eeg{1,trialno}(channel,:);

        plot(singleeeg)
        %ylim([min(singleeeg) max(singleeeg)]);
         ylim([-400 400]);
         xlim([0 length(singleeeg)]);
        ylabel(sprintf('Tetrode %s',num2str(channel),num2str(trialno)));
        set(gca,'xtick',[],'ytick',[]);
    end
end
end





% end
%powerspectrum=[powerspectrum;S2];
%hold on;



%semilogx(freqs,S2);
%end




% initiate parameters
%
%
%     animalname = out.animal{1,1};                              % to label graphs later
%     eventtype = 'eeg';                              % must be the same name as the events' data structure
%     days = max(out.epochs{1,1}(:,1));
%
%     ref_tetrode = 1;                                   % reference tetrode
%     tetrodes = channel; %ceil(channel/4);                                    % tetrodes to analyze
%     tetno=max(tetrodes);
%     epno=max(out.epochs{1,1}(:,2));
%
%
% for ii= 1; %1:size(out.epochs{1,1},1)
%
%
%     d=days;% days to analyze
%     t=tetrodes;
%
%     % epochs to analyze
%     % velocity_state_1 = 2;                               % state 1 (immobile)
%     % velocity_state_2 = 8;                               % state 2 (run)
%     windowsize_sec = 20;                               %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            % window length
%     Fs = 1500;                                          % sample rate of EEG
%     halfwindow_samp = (windowsize_sec*1500)/2;
%
%
%     filename='';
%
%     epochs = out.epochs{1,1}(ii,2);
%     e=epochs;
%     eventno=length(outdata{1,ii}.eeg);
%
%     % make cells
% %     eegwindows=cell(days(end),epno,tetno,eventno);
% %     allwindowedspectrograms=cell(days(end),epno,tetno,eventno);
% %     continuousspectrograms=cell(days(end),epno,tetno,eventno);
% %     meandayspectra=cell(days(end),tetno,eventno);
% %     stddayspectra=cell(days(end),tetno,eventno);
% %     zscorespectrograms=cell(days(end),epno,tetno,eventno);
% %     meanspectrograms_epoch=cell(days(end),epno,tetno,eventno);
%
%
%     eegwindows={};
%     allwindowedspectrograms={};
%     continuousspectrograms={};
%     meandayspectra={};
%     stddayspectra={};
%     zscorespectrograms={};
%     meanspectrograms_epoch={};
%
%
%
%         %% Set these parameters manually.
%         % Chronux params are specified below (fourth block)
%
%         %% First, collects all event timestamps from chosen ref. tetrode.
%
%
%         %% Second, copy the raw eeg data into windows aligned to the event starttime.
%
%         % load eeg struct
%
%         eegstruct=loadeegstruct(out.animal{1,2},animalname,'eeg',days,epochs,tetrodes);
%
%
%
%
%
%         %% Fourth, calculate all individual spectrograms.
%
%         % holds all events' spectograms
%
%         % chronux params
%         movingwin = [1000 100]/1000;
%         params.Fs = 1500;
%         params.err = [2 0.05];
%         params.fpass = [0 25];
%         params.tapers = [3 5];
%
%
%          for en= 1:size(outdata{1,ii}.eeg,2)
%
%         eegwindows{d}{e}{t}{en}=outdata{1,ii}.eeg{1,en}(t,:);
%          end
%
%
%
%
%         for en=1:eventno
%             flag=0;
%             %for r=1:size(eegwindows{d}{e}{t}{en},1)
%                 [S,times,freqs,Serr]=mtspecgramc(eegwindows{d}{e}{t}{en},movingwin,params);
% %                 if flag==0                              % this if clause initializes the 3D matrix of all spectrograms for a given tetrode
% %                     allwindowedspectrograms{d,e,t,en}=nan(size(S,1),size(S,2),size(eegwindows{d,e,t,en},1));
% %                     flag=1;
% %                 end
%                 allwindowedspectrograms{d}{e}{t}{en}=S;   % adds the S variable (2D matrix) to the third matrix dimension
%             %end
%         end
%
%
%         %% Fifth, calculates the continuous spectrogram, day mean, and std spectra for each tetrode.  (long calculation)
%
%
%
%         dummy=[];
%
%         [S_full,junkt,junkf,junkserr] = mtspecgramc(eegstruct{d}{e}{t}.data,movingwin,params);
%         dummy=[dummy;S_full];
%         continuousspectrograms{d}{e}{t}{en}=S_full;
%
%         meandayspectra{t}=mean(dummy,1);
%         stddayspectra{t}=std(dummy,1);
%
%         % z-score the continuous spectrogram
%
%         for en=eventno
%             continuousspectrograms{d}{e}{t}{en}=bsxfun(@minus,continuousspectrograms{d}{e}{t}{en},meandayspectra{t});
%             continuousspectrograms{d}{e}{t}{en}=bsxfun(@rdivide,continuousspectrograms{d}{e}{t}{en},stddayspectra{t});
%         end
%
%         disp('.finished with computing spectrograms.')
%         %save(filename,'-7.3');
%
%         %% Sixth, z-scores the individual spectrograms ('zscorespectrograms').
%
%         for en=1:eventno;
%             if ~isempty(allwindowedspectrograms{d}{e}{t}{en})       % in case there are no events
%                 for r=1:size(allwindowedspectrograms{d}{e}{t}{en},3)
%                     zscorespectrograms{d}{e}{t}{en}=bsxfun(@minus,allwindowedspectrograms{d}{e}{t}{en},meandayspectra{t});
%                     zscorespectrograms{d}{e}{t}{en}=bsxfun(@rdivide,allwindowedspectrograms{d}{e}{t}{en},stddayspectra{t});
%                 end
%             end
%         end
%
%
%
% %         %% Seventh, calculates mean spectrograms: (1) state, (2) epochs (3) overall.
% %
% %         %meanspectrograms_state=cell(days(end),tetno,3);             %% discards epoch, pools all within state
% %         % discards states, pools all within epoch
% %         % meanspectrograms_allevents=cell(days(end),tetno);           %% pools all events
% %
% %         %epochs=[2];
% %
% %         for en=1:eventno;
% %         dummy=[];
% %
% %             if ~isempty(zscorespectrograms{d}{e}{t}{en})
% %                 dummy=cat(3,dummy,zscorespectrograms{d}{e}{t}{en});
% %             end
% %
% %         %meanspectrograms_state{d,t,s}=mean(dummy,3);
% %         %meanspectrograms_allevents{d,t,en}=cat(3,dummy,meanspectrograms_allevents{d,t,en});
% %
% %         %meanspectrograms_allevents{d,t,en}=mean(meanspectrograms_allevents{d,t,en},3);
% %
% %         dummy=[];
% %
% %         if ~isempty(zscorespectrograms{d}{e}{t}{en})
% %             dummy=cat(3,dummy,zscorespectrograms{d}{e}{t}{en});
% %         end
% %
% %         meanspectrograms_epoch{d}{e}{t}{en}=mean(dummy,3);
% %         end
%
%
%
%         %% (ALTERNATIVE) Plot mean spectrogram for each epoch.
%
%         %         for d=days
%         %             for e=epochs
%         %                 for t=tetrodes
%         %                     %h = subplot(1,11,t);
%         %                     figure;
%         %                     imagesc(times,freqs,meanspectrograms_epoch{d,e,t,en}',[-0.2,5]);
%         %                     set(gca,'YDir','normal');
%         %                     string = sprintf('%d',t);
%         %                     %title(string);
%         %                     colormap('jet')
%         %                     colorbar
%         %
%         %                     %%% TITLE
%         %                     % eventcount = sum(statecounts(:,e));
%         %                     h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
%         %                         1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%         %                     text(0.5, 1,sprintf('Mean %s spectrogram, %s, Day %d, reftet=%d, epoch %d ', ...
%         %                         eventtype,animalname,d,ref_tetrode,e),'HorizontalAlignment'...
%         %                         ,'center','VerticalAlignment', 'top');
%         %
%         %                     %     % animal_day_epoch
%         %                     cd(strcat(f.animal{1,2},'Plot/'));
%         %                     figurename = strcat(animals{1,1},'_eegspec_reward',num2str(d),'_e',num2str(e),'_t',num2str(t),'_w',num2str(jj));
%         %                     %
%         %                     %saveas(gcf, figurename, 'pdf');
%         %                     %
%         %                     %Closes the figure
%         %                     close;
%         %                     %
%         %                 end
%         %             end
%         %         end
%     end
%
%
% % Tenth, plot some random individual event spectrograms.
%
% % set these manually to your taste
% d=days;
% e=2;                                 % epoch manually
% events_plot=10;                         % no. events to plot
% tetrodes_plot=channel;                  % tetrodes to plot
% %states_plot=1;
% %
%
% t=tetrodes;
% %for i=1:events_plot
% % N=size(zscorespectrograms{d,e,t,en},3);        % # events to choose from
% % eventno=ceil(N*rand);
% figure
% %title(string);
%
%
% for i=1:events_plot
%     N=size(zscorespectrograms{d}{e}{t},2);       % # events to choose from
%     en=ceil(N*rand);
%     %eventno=2;
%
%     for t=tetrodes_plot
%         figure;
%         %subplot(1,8,t-7)
%         imagesc(times-windowsize_sec/2,freqs,zscorespectrograms{d}{e}{t}{en}',[-0.5 5]);
%         colormap('jet')
%         set(gca,'YDir','normal');
%         %string = sprintf('%d',t);
%         %title(string);
%     end
% end
% %%% TITLE code.
% h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
%     1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 1,sprintf('Single %s event spectrogram, %s, Day %d, epoch %d reftet=%d (event #%d)', ...
%     eventtype,animalname,d, e,ref_tetrode,eventno),'HorizontalAlignment'...
%     ,'center','VerticalAlignment', 'top');
% %end.
% cd(strcat(f.animal{1,2},'Plot/'));
% figurename = strcat(animals{1,1},'_eegsingle_reward',num2str(d),'_e',num2str(e),'_t',num2str(t));
% %
% %saveas(gcf, figurename, 'pdf');
% %
% %Closes the figure
% close;

% plot unique cells

% compare firing rate to trajectory distance



%     % ----Saving----
%     % Saves figure as pdf
%     % First checks if a folder called Plot exists in the processed data folder,
%     % if not, creates it.
%
%     cd(f.animal{1,2});
%     plotdir = dir('Plot');
%     if (isempty(plotdir))
%         %an a plot folder needs to be created
%         !mkdir Plot
%     end
%
%     % change to that directory and saves the figure with file name
%     % animal_day_epoch
%     cd(strcat(f.animal{1,2},'Plot/'));
%     figurename = strcat(animals{1,1},'_intertrialspikelocation_d',num2str(daytext),'_e',num2str(dayindex),'_t',num2str(uniquecells(ii,1)),'_c',num2str(uniquecells(ii,2)));
%
%     saveas(gcf, figurename, 'pdf');
%
%     %Closes the figure
%     close;
%
%
%
%
%     end
% end
%end
