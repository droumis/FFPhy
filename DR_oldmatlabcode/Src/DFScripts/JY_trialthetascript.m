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

epochfilter{1} = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch,  4)'];

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

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);
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

out=setfilterfunction(f, 'JY_trialtheta', {'data','theta'});

out=runfilter(out);

outdata=out.output{1,1};


% generate index to show day epoch and tetrode channel



signalstd=std(double(out.output{1,1}{1,1}.intertrialthetaphase{1,1}(15:21,:)),[],1);

plot(signalstd);

data1=double(out.output{1,1}{1,1}.thetaphase{1,2}(2,:)');
data2=double(out.output{1,1}{1,1}.thetaphase{1,2}(3,:)');
params = [];
params.Fs = 1500;
params.fpass = [1 25];
params.err = [1 0.05];
params.tapers = [3 5];
window =[1000 10]/1000;

[C,phi,S12,S1,S2,t,f]=cohgramc(data1,data2,window,params);


 figure;
                    imagesc(t,f,C');
                    set(gca,'YDir','normal');
                    string = sprintf('%d',times);
                    %title(string);
                    colormap('jet')
                    colorbar
                    

figure;
channel=2;
% calculate power for all epochs and average
powerspectrum=[];
for trialno=1:size(out.output{1,1}{1,1}.theta,2)

% for i=1:size(out.output{1,1}{1,1}.eeg{1,trialno},1)
%     h=subplot(size(out.output{1,1}{1,1}.eeg{1,trialno},1),1,i);
%     plot(out.output{1,1}{1,1}.eeg{1,trialno}(i,:));
%     set(h,'XTick',[],'YTick',[]);
%     ylabel(sprintf('%s',num2str(i)));
% end



singleeeg=out.output{1,1}{1,1}.eeg{1,trialno}(channel,:);

 movingwin = [1000 50]/1000;
        params.Fs = 1500;
        params.err = [2 0.05];
        params.fpass = [0 200];
        params.tapers = [3 5];


        %[S,times,freqs,Serr] = mtspecgramc(singleeeg,movingwin,params);

     %figure;

        %imagesc(times,freqs,S',[-0.2 max(S(:))]);

        %imagesc(times,freqs,S',[-0.2 100]);

   %set(gca,'YDir','normal');

   %figure;

   [S2,freqs,Serr] = mtspectrumc(singleeeg,params);
   
   %powerspectrum=[powerspectrum;S2];
   hold on;
   
   
   semilogx(freqs,S2);
end




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
