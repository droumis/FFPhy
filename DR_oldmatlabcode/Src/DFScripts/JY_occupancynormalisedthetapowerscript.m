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

%epochfilter{1} = [''];

%epochtype='Run';

%epochfilter{1} = ['isequal($epochtype, ''Run'')'];
epochfilter = ['isequal($epoch, 6 )'];

%cellfilter = '(isequal($area, ''ACC'') && ($meanrate >0 ) )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0 ) )'  ;
timefilter = {{ 'JY_getlinvelocity', '$velocity  <3' } };
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

%timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%    {'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))},{'JY_getbarrier','($barrier== 0)'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { 'JY_getlinvelocity', '$velocity > 0' };


tetrodefilter=['isequal($area, ''HP'')'];

f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'excludetimefilter', timefilter,'eegtetrodes',tetrodefilter);

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
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

out=setfilterfunction(f, 'JY_occupancynormalisedthetapower', {'data','meaneegspectrograms','theta',...
    'tetinfo'});


out=runfilter(out);

outdata=out.output{1,1};


%% plotting

for ii=1:size(out.epochs{1,1},1)
    
    animalname = out.animal{1,1};                              % to label graphs later
                           % must be the same name as the events' data structure
        days = out.epochs{1,1}(ii,1);                                           % days to analyze
                             % tetrodes to analyze
        
        
        dataindex=cell2mat(f.eegdata{1,1}');
          
        
        for iii=1:size(dataindex,1);  
        
            
        epoch = dataindex(iii,1);
        tetrode= dataindex(iii,2);
        
        data=outdata(1,iii).epochmeanzscore;
        times=outdata(1,iii).times;
        freq=outdata(1,iii).freq;
        
                    %h = subplot(1,11,t);
                    figure;
                    imagesc(times,freq,data');
                    set(gca,'YDir','normal');
                    string = sprintf('%d',times);
                    %title(string);
                    colormap('jet')
                    colorbar
                    
                    %%% TITLE
                    % eventcount = sum(statecounts(:,e));
                    h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
                        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
                    text(0.5, 1,sprintf('Mean reward triggered spectrogram, %s, day %d, epoch %d, tetrode %s', ...
                        animalname,days,epoch,num2str(tetrode)),'HorizontalAlignment'...
                        ,'center','VerticalAlignment', 'top','FontSize',8);
                    
                    %     % animal_day_epoch
                    cd(strcat(f.animal{1,2},'Plot/'));
                    figurename = sprintf('%seegspecrwd%s_%s_%s',animals{1,1},num2str(days),num2str(epoch),num2str(tetrode));
                    %
                    
                    print(figurename,'-depsc');
                    
                    %saveas(gcf, figurename, 'pdf');
                    %
                    %Closes the figure
                    close;
                    %
                end
            end