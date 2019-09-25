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

days='[5]';%,'1:10';
%days = '[1:1]';
%days = '[9:9]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = ['1:7'];

%epochtype='Run';
%epochfilter = ['isequal($epochtype, ''Sleep'')'];
%epochfilter = ['isequal($epochtype, ''Run'')'];
epochfilter = ['isequal($epoch,  6)'];
%epochfilter = ['ismember($epoch,  1:7)'];

%cellfilter = '(isequal($area, ''ACC'') && ($meanrate >0 ) )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0 ) )'  ;

%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

%timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'}};

%exclude all times with velocity
timefilter = {{ 'JY_getlinvelocity', '$velocity  >3' } };

%timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%    {'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))},{'JY_getbarrier','($barrier== 0)'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { 'JY_getlinvelocity', '$velocity > 0' };

tetrodefilter=['isequal($area, ''HP'')'];
tetrodereference=[18];

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter, 'excludetimefilter', timefilter,'eegtetrodes',tetrodefilter);
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'excludetimefilter', timefilter,'eegtetrodes',tetrodefilter);
%f = JY_createfilter(days,'animal',animals,'days',days,'excludetimefilter', timefilter,'eegtetrodes',tetrodefilter);

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);
%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'JY_eeganal';

f = setfilteriterator(f,iterator);

%out=setfilterfunction(f, 'JY_alignspectrumtoripple',{'eeg','meaneegspectrograms'});
%out=setfilterfunction(f, 'JY_alignspectrumtoripple',{'eeg','meaneegspectrograms','ripples','tetinfo'},'tetrodereference',[17]);
out=setfilterfunction(f, 'JY_calcriptriggeredspectrogram',{'data','eeg','meaneegspectrograms',...
    'ripples','ripple','tetinfo'},'tetrodereference',[18]);
out=runfilter(out);

outdata=out.output{1,1};

% process each epoch separately

% based on event_spectrograms.m from KKay


% generate index to show day epoch and tetrode channel
% % 
for ii=1:size(out.epochs{1,1},1)
    
    animalname = out.animal{1,1};                              % to label graphs later
    % must be the same name as the events' data structure
    days = out.epochs{1,1}(ii,1);                                           % days to analyze
    % tetrodes to analyze
    
    
    %dataindex=cell2mat(f.eegdata{1,1}');
    dataindex=f.eegdata{1,1}{1,ii};
    
    for iii=1:size(dataindex,1);
        
        epoch = dataindex(iii,1);
        tetrode= dataindex(iii,2);
        
        data=outdata(1,iii).epochmeanzscore;
        times=outdata(1,iii).times;
        freq=outdata(1,iii).freq;
        
        % plot the ripple eeg
        
        meanrippleeeg=outdata(1,iii).meanrippleeeg;
        
        figure; set(gcf,'position',[0 0 600 600]); set(gcf,'PaperPositionMode','auto');
        
        nrows=2;
        columns=1;
        %gap_h=0.0000000001;
        %gap_w=0.00000001;
        gap_h=0.005; gap_w=0.05;
        marg_h=[0.05 0.10]; marg_w=[0.1 0.1];
        
        ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
        axes(ha(1));
        
        plot(meanrippleeeg);
        xlim([0 length(meanrippleeeg)]); title(sprintf('Mean ripple'));
        colorbar;
        
        %h = subplot(1,11,t);
        axes(ha(2))
        imagesc(times,freq,data',[-0.5 0.5]);
        set(gca,'YDir','normal');
        string = sprintf('%d',times);
        %title(string);
        colormap('jet')
        colorbar;
        
        %%% TITLE
        % eventcount = sum(statecounts(:,e));
        h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
            1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        %text(0.5, 1,sprintf('Mean %s spectrogram, %s, Day %d, reftet=%d, epoch %d ', ...
        %    eventtype,animalname,d,ref_tetrode,e),'HorizontalAlignment'...
        %    ,'center','VerticalAlignment', 'top');
              text(0.5, 1,sprintf('Mean ripple triggered spectrogram, %s, day %d, epoch %d, tetrode %s', ...
                        animalname,days,epoch,num2str(tetrode)),'HorizontalAlignment'...
                        ,'center','VerticalAlignment', 'top','FontSize',8);

        %     % animal_day_epoch
        cd(strcat(f.animal{1,2},'Plot/'));
        %figurename = strcat(animals{1,1},'_eegspec_reward',num2str(d),'_e',num2str(e),'_t',num2str(t),'_w',num2str(jj));
        %
           figurename = sprintf('%seegspecrpl%s_%s_%s',animals{1,1},num2str(days),num2str(epoch),num2str(tetrode));
                    %
                    
                    %print(figurename,'-depsc');

        %saveas(gcf, figurename, 'pdf');
        %
        %Closes the figure
        %close;
        %
    end
end


%% pool epoch results for each tetrode

% get indices

% epochindices=cell2mat(cellfun(@(x) x(:,2)',out.eegdata{1,1},'UniformOutput',false))';
% 
% uniquetet=unique(epochindices);
% 
% for ii=uniquetet'
%     
%     eegindices=find(epochindices==ii);
%     
%     meanrawz=arrayfun(@(x) x.epochrawzscore,outdata(eegindices),'UniformOutput',false);
%     meanrawz=cat(3,meanrawz{:});
%     meanz=mean(meanrawz,3);
%     
%     times=outdata(1).times;
%     freq=outdata(1).freq;
%     
%     figure;
%     imagesc(times,freq,meanz');
%     set(gca,'YDir','normal');
%     string = sprintf('%d',times);
%     %title(string);
%     colormap('jet')
%     colorbar
% end


% % Tenth, plot some random individual event spectrograms.
%
% % set these manually to your taste
% d=days;
% e=epochs;                                 % epoch manually
% events_plot=10;                         % no. events to plot
% tetrodes_plot=19;                  % tetrodes to plot
% %states_plot=1;
% %
%
% t=19;
% %for i=1:events_plot
% % N=size(zscorespectrograms{d,e,t},3);        % # events to choose from
% % eventno=ceil(N*rand);
% figure
% %title(string);
% for i=1:events_plot
%      N=size(zscorespectrograms{d,e,t},3);        % # events to choose from
%             eventno=ceil(N*rand);
%
% for t=tetrodes_plot
%     figure;
%     %subplot(1,8,t-7)
%     imagesc(times-windowsize_sec/2,freqs,zscorespectrograms{d,e,t}(:,:,eventno)',[-0.2,4]);
%     colormap('jet')
%     set(gca,'YDir','normal');
%     %string = sprintf('%d',t);
%     %title(string);
% end
% end
% %%% TITLE code.
% h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
%     1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 1,sprintf('Single %s event spectrogram, %s, Day %d, reftet=%d (event #%d)', ...
%     eventtype,animalname,d,ref_tetrode,eventno),'HorizontalAlignment'...
%     ,'center','VerticalAlignment', 'top');
% %end.
% cd(strcat(f.animal{1,2},'Plot/'));
% figurename = strcat(animals{1,1},'_eegsingle_reward',num2str(d),'_e',num2str(e),'_t',num2str(t));
% %
% saveas(gcf, figurename, 'pdf');
% %
% %Closes the figure
% close;
%
% % plot unique cells
%
% % compare firing rate to trajectory distance
%
%
%
% %     % ----Saving----
% %     % Saves figure as pdf
% %     % First checks if a folder called Plot exists in the processed data folder,
% %     % if not, creates it.
% %
% %     cd(f.animal{1,2});
% %     plotdir = dir('Plot');
% %     if (isempty(plotdir))
% %         %an a plot folder needs to be created
% %         !mkdir Plot
% %     end
% %
% %     % change to that directory and saves the figure with file name
% %     % animal_day_epoch
% %     cd(strcat(f.animal{1,2},'Plot/'));
% %     figurename = strcat(animals{1,1},'_intertrialspikelocation_d',num2str(daytext),'_e',num2str(dayindex),'_t',num2str(uniquecells(ii,1)),'_c',num2str(uniquecells(ii,2)));
% %
% %     saveas(gcf, figurename, 'pdf');
% %
% %     %Closes the figure
% %     close;
% %
% %
% %
% %
% %     end
% % end
% %end
