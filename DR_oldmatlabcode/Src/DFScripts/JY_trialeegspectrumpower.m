%plots the spectral power, averaged over all trials.




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
animals = {'CML21'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filterionno

days='[16]';%,'1:10';
%days = '[1:1]';
%days = '[9:9]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = [''];

epochtype='Run';

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

out=setfilterfunction(f, 'JY_trialeeg', {'data','eeg'});

out=runfilter(out);

outdata=out.output{1,1};

% process each epoch separately

% based on event_spectrograms.m from KKay







for jj=1:size(out.epochs{1,1},1)
     figure;
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[2500 1000]);
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'position',[10 10 2300 800]);


    day=out.epochs{1,1}(jj,1);
    epoch=out.epochs{1,1}(jj,2);
    
    % generate index to show day epoch and tetrode channel
    plotchannels=[1 2 4 5 6 7 8 9 10 11 12 13 19 20 21];
    noofplots=length(plotchannels);
    movingwin = [1000 100]/1000;
    params.Fs = 1500;
    params.err = [2 0.05];
    params.fpass = [0 400];
    params.tapers = [3 5];
    params.trialave=1;
    
    nrows=4;
    columns=4;
    gap_h=0.01;
    gap_w=0.01;
    marg_h=[0.1 0.1];
    marg_w=[0.01 0.01];
    ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
    ax=1;
    
   
    for ii=1:length(plotchannels);
        
        channel=plotchannels(ii);

        chunkeddata=[];
        chunkedinterdata=[];
        
        % calculate power for all epochs and average
        powerspectrum=[];

        for trialno=1:size(out.output{1,1}{1,jj}.eeg,2)
            
            singleeeg=out.output{1,1}{1,jj}.eeg{1,trialno}(channel,:);
            
            % break data into 1s chunks
            
            maxsize=floor(size(singleeeg,2)/params.Fs);
            tempchunkeddata=reshape(singleeeg(1:maxsize*params.Fs),params.Fs,maxsize);
            
            
            chunkeddata=[chunkeddata tempchunkeddata];
        end
        
        % calculate intertrial power for all epochs and average
        powerspectrum=[];

        for intertrialno=1:size(out.output{1,1}{1,jj}.intertrialeeg,2)
            
            singleintertrialeeg=out.output{1,1}{1,jj}.intertrialeeg{1,intertrialno}(channel,:);
            
            % break data into 1s chunks
            
            maxsize=floor(size(singleintertrialeeg,2)/params.Fs);
            tempchunkeddata=reshape(singleintertrialeeg(1:maxsize*params.Fs),params.Fs,maxsize);
            
            
            chunkedinterdata=[chunkedinterdata tempchunkeddata];

        end     

        %for i=1:size(out.output{1,1}{1,jj}.eeg{1,trialno},1)
        %         h=subplot(size(out.output{1,1}{1,jj}.eeg{1,trialno},1),1,i);
        %         plot(out.output{1,1}{1,jj}.eeg{1,trialno}(i,:));
        %         set(h,'XTick',[],'YTick',[]);
        %         ylabel(sprintf('%s',num2str(i)));
        %     end
        
        [S2,freqs,Serr] = mtspectrumc(chunkeddata,params);
        [S2i,freqsi,Serri] = mtspectrumc(chunkedinterdata,params);
        
        axes(ha(ax));
        
        %semilogx(freqs,S2);
        hold on;
        
%         semilogx(freqs,Serr(1,:),'.r');
%         semilogx(freqs,Serr(2,:),'.r');

        %semilogx(freqsi,S2i,'g');
        %hold on;
        
%         semilogx(freqsi,Serri(1,:),'.c');
%         semilogx(freqsi,Serri(2,:),'.c');


        %figure;
        bar(freqsi,S2./S2i);
        set(gca,'YScale','log')


        xlim(params.fpass);
        set(gca,'ytick',[]);
        set(gca,'YScale','log')
   
        % print text
        xlimval=xlim;
        ylimval=ylim;
        % distlabel=[0.35*(xlimval(2)-xlimval(1))+xlimval(1) 0.05*(ylimval(2)-ylimval(1))+ylimval(1)];
        % text(distlabel(1),distlabel(2),sprintf('day %s epoch %s',num2str(d),num2str(e)),'FontSize',12,'Color','k');
        % print number of ripples
        distlabel2=[0.10*(xlimval(2)-xlimval(1))+xlimval(1) 0.85*(ylimval(2)-ylimval(1))+ylimval(1)];
        text(distlabel2(1),distlabel2(2),sprintf('t%s',num2str(channel)),'FontSize',12,'Color','k');
        ax=ax+1;
        
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    
    
    text(0.5, 0.98,sprintf('PSD for %s Day %s Epoch %s', animals{1,1},...
        num2str(day), num2str(epoch)),'FontSize',12,'HorizontalAlignment',...
        'center','VerticalAlignment', 'top');
    
    % ----Saving----
    % Saves figure as pdf
    % First checks if a folder called Plot exists in the processed data folder,
    % if not, creates it.
      set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[30 15]);

    cd(f.animal{1,2});
    plotdir = dir('Plot');
    if (isempty(plotdir))
        %an a plot folder needs to be created
        !mkdir Plot
    end
    
    % change to that directory and saves the figure with file name
    % animal_day_epoch
    cd(strcat(f.animal{1,2},'Plot/'));
    figurename = strcat(animals{1,1},'_spectrumpower_d',num2str(day),'_e',num2str(epoch));
    
     saveas(gcf, figurename, 'fig');
     saveas(gcf, figurename, 'pdf');
    
    %Closes the figure
    close;
    
end
