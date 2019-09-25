
% plots all coherence pairs between tetrodes, with each epoch as a single
% sample

runscript = 0
% pre parameters
animal_list = {'Frank'}
windowsize = 0.5;  % in seconds
    % region 1
area1 = 'CA1';
numcells1 = 2;
    tetfilter1 = ['(isequal($area, ''' area1 ''') && ($numcells >= ''' numcells1 '''))'];
    % region 2
area2 = 'CA3';
numcells2 = 2;
    tetfilter2 = ['(isequal($area, ''' area2 ''') && ($numcells >= ''' numcells2 '''))'];

epoch = 'run';
timefilt_flag = 'all';  

% post parameters
freqrange = [70 90];


if runscript

% Animal Selection
animals = animal_list;

% dayfilter
%dayfilter = 3:8; 

% Epoch Filter
if strcmp(epoch,'run')
    epochfilter = '(isequal($type, ''run''))';
elseif strcmp(epoch,'sleep')
    epochfilter = '(isequal($type, ''sleep''))';
elseif strcmp(epoch,'all')
    epochfilter = '(isequal($type, ''run'') || isequal($type, ''sleep''))';
else
    error('specify epoch type!')
end

% Time Filter
if strcmp(timefilt_flag,'all')
    timefilter = { {'kk_get2dstate', '($velocity >= 0)'} };
elseif strcmp(timefilt_flag,'theta')
    timefilter = { {'kk_getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
                   {'kk_get2dstate', '($velocity >= 4)'} };
elseif strcmp(timefilt_flag,'rip')
    timefilter = { {'kk_getriptimes', '($nripples == 2)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...
                    {'kk_get2dstate', '($velocity < 4)'} };                
end

% Iterator
iterator = 'kk_singleepochanal';

% Filter Creation
%f = createfilter('animal', animals, 'epochs', epochfilter, 'days',dayfilter, 'excludetime', timefilter, 'iterator', iterator);
f = createfilter('animal', animals, 'epochs', epochfilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
f = setfilterfunction(f, 'dfakk_getcoherencepairs', {'eeg','tetinfo'},freqrange,tetfilter1,tetfilter2,'windowsize',windowsize);

% Run Analysis
f = runfilter(f);

end


% Parameter string (for figures).
%params_string = [epoch ' epochs / ' 'timefilter: ' timefilt_flag ' / minthresh: ' num2str(minthresh) ...
%                 ' / ngammal: ' num2str(ngammal) ' / phaseregion: ' phasereg ' / onephasetet: ' num2str(onephasetet)]
                            

% Plot all coherence values in a histogram

pairs = [];
coherograms = [];
freqflag = 0;

for ep = 1:length(f.output{1})
   
    pairs = [pairs ; f.output{1}(ep).pairs];
    coherograms = [coherograms ; f.output{1}(ep).coherograms];
    if ~isempty(f.output{1}(ep).freqs)
        freqs = f.output{1}(ep).freqs;
        lowindex = lookup(freqrange(1),freqs);
        highindex = lookup(freqrange(2),freqs);
    end
    
end

coherences = [];

for c = 1:size(coherograms,1)
 
    coherences(c) = mean(coherograms(c,lowindex:highindex));
 
end


% plot etc.

figure
plot(f.output{1}(1).freqs,coherograms)   % raw coherograms
figure
plot(f.output{1}(1).freqs,mean(coherograms,1))  % grand mean coherogram
figure
hist(coherences,200)   % histogram


% example coherograms
%for i = 150:200
%    figure
%    plot(f.output{1}(20).freqs,coherograms(i,:))
%    ylim([0 1])
%end

