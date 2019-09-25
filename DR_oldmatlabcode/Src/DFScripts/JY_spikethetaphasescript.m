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

days='[14]';%,'1:10';
%days = '[1:1]';
%days = '[9:9]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = [''];

epochtype='Run';

epochfilter{1} = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch,  4)'];

cellfilter = '(isequal($area, ''ACC'') && ($meanrate >0 ) )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
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


%get filtered theta and eeg
%--------------------------------------------
iterator = 'JY_singleepochanal';

f = setfilteriterator(f,iterator);

out=setfilterfunction(f, 'JY_trialtheta', {'data','theta'});

out=runfilter(out);

outdata=out.output{1,1};



eeg=setfilterfunction(f, 'JY_trialeeg', {'data','eeg'});

eeg=runfilter(eeg);

outeeg=eeg.output{1,1};



spikelocation = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter);
%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';

spikelocation = setfilteriterator(spikelocation,iterator);

outspike=setfilterfunction(spikelocation, 'JY_trialspikelocation', {'spikes','data','linpos'});

outspike=runfilter(outspike);


%plotting

thetatetrode=13;

epochindex=2;

trialno=6;

%plot filtered theta with eeg

eegtime=[outeeg{1,epochindex}.trialtimes(trialno,1):10000/outeeg{1,epochindex}.samplerate:outeeg{1,epochindex}.trialtimes(trialno,2)];

eeg=outeeg{1,epochindex}.eeg{1,trialno}(thetatetrode,:);

thetatime=[outdata{1,epochindex}.trialtimes(trialno,1):10000/outdata{1,epochindex}.samplerate:outdata{1,epochindex}.trialtimes(trialno,2)];

theta=outdata{1,epochindex}.theta{1,trialno}(thetatetrode,:);


%plot(eegtime,eeg(1:size(eegtime,2)));

hold on;

plot(thetatime,theta(1:size(thetatime,2)),'k');


% put spikes on top

spikes=outspike.output{1,1}{1,epochindex}.trajectory{1,trialno}.spikes(:,1)*10000;


hold on;

thetaspike=lookup(spikes,thetatime);


plot(thetatime(thetaspike),theta(thetaspike),'.r');


% look up theta phase
allspikes=outspike.output{1,1}{1,epochindex}.allspiketimes;
thetaphase=outdata{1,epochindex}.allthetaphase(thetatetrode,:);
thetaphasetime=outdata{1,epochindex}.timeindex(thetatetrode,:);
spikethetaphase=double(thetaphase(lookup(allspikes,thetaphasetime)));

hist(spikethetaphase/10000,20);













