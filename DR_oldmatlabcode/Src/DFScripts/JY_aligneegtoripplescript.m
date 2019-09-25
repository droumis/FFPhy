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

%epochfilter = ['isequal($epochtype, ''Run'')'];
epochfilter = ['isequal($epoch,  4)'];

%cellfilter = '(isequal($area, ''ACC'') && ($meanrate >0 ) )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0 ) )'  ;

%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

%timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'}};

%exclude all times with velocity
%timefilter = {{ 'JY_getlinvelocity', '$velocity  >3' } };

%timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%    {'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))},{'JY_getbarrier','($barrier== 0)'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { 'JY_getlinvelocity', '$velocity > 0' };

tetrodefilter=['isequal($area, ''ACC'')'];
tetrodereference=[16];

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter, 'excludetimefilter', timefilter,'eegtetrodes',tetrodefilter);
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'eegtetrodes',tetrodefilter);

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);
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

%out=setfilterfunction(f, 'JY_alignspectrumtoripple',{'eeg','meaneegspectrograms'});
out=setfilterfunction(f, 'JY_aligneegtoripple',{'data','eeg','meaneegspectrograms','ripple','ripples'},'tetrodereference',[18],'std',[5]);
%out=setfilterfunction(f, 'JY_calcriptriggeredspectrogram',{'eeg','meaneegspectrograms','ripples','tetinfo'},'tetrodereference',[17]);
out=runfilter(out);

outdata=out.output{1,1};

% process each epoch separately

% based on event_spectrograms.m from KKay


% generate index to show day epoch and tetrode channel
%% plot raw eeg trace trial

% generate index to show day epoch and tetrode channel
for trialno=1:size(out.output{1,1}{1,1}.ripplespeed,2);
%figure;
tetrodes=[1:21];
speed=out.output{1,1}{1,1}.ripplespeed{trialno};
ripple=out.output{1,1}{1,1}.rippleripple{trialno};
ripplesize=out.output{1,1}{1,1}.ripplesize(trialno);
figure; set(gcf,'position',[500 200 300 800]); set(gcf,'PaperPositionMode','auto');

nrows=size(tetrodes,2)+1;
columns=1;
%gap_h=0.0000000001;
%gap_w=0.00000001;
gap_h=0.005; gap_w=0.05;
marg_h=[0.05 0.10]; marg_w=[0.1 0.1];

ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
axes(ha(1));
plot(speed);
xlim([0 length(speed)]); title(sprintf('Trial %s %s std',num2str(trialno), num2str(ripplesize)));

axes(ha(2));
plot(ripple)
xlim([0 length(ripple)]);
%singleeeg=out.output{1,1}{1,1}.intertrialeeg{1,trialno}(channel,:);
%speed=out.output{1,1}{1,1}.intertrialspeed{trialno};

for ii=3:nrows;
    for channel=tetrodes(ii-1);
        ax=ii;
        axes(ha(ax));
        
        singleeeg=out.output{1,1}{1,1}.rippleeeg{1,trialno}(channel,:);

        plot(singleeeg)
        %ylim([min(singleeeg) max(singleeeg)]);
         ylim([-300 300]);
         xlim([0 length(singleeeg)]);
        ylabel(sprintf('Tetrode %s',num2str(channel),num2str(trialno)));
        set(gca,'xtick',[],'ytick',[]);
    end
end
end
