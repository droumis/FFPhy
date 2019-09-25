%Animal selection
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = '($dailyexposure == 1) & isequal($description, ''TrackB'')';
epochfilter{2} = '($dailyexposure == 1) & isequal($description, ''TrackA'')';

cellfilter = '($meanrate<7) & (isequal($area,''CA1'') | isequal($area,''CA3''))';

%Define iterator
iterator = 'multicellanal';

%% RUN TRAINING FILTER

%Define timefilter for training data
timefilter = { {'getriptimes','($nripples == 0)', [], 'cellfilter', ...
'((isequal($area, ''CA3''))|(isequal($area, ''CA1'')))'} };

%create training data by calulating the linearized rates of all cells
trainingfilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);

%Options set: bin size = 2 cm, peakthresh = 3 std above the mean
trainingfilter = setfilterfunction(trainingfilter, 'calcpopulationlinfields', {'spikes','linpos'},2,3);
trainingfilter = runfilter(trainingfilter);

%% RUN DECODING FILTER FOR RUNS

%Filter creation for position decoding
epochfilter = [];
epochfilter{1} = 'isequal($description, ''TrackB'')';
epochfilter{2} = 'isequal($description, ''TrackA'')';

%Define minimum number of cells for candidate replay events
cellcountthresh = 5;

timefilter = {{'get2dstate', '((abs($velocity) < 4))'}};
clear decodefilter
global decodefilter ;
decodefilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
decodefilter = setfilterfunction(decodefilter, 'getpopulationevents', {'spikes','linpos', 'pos','ripples','cellinfo'},cellcountthresh);
decodefilter = runfilter(decodefilter);

%% RUN DECODING FILTER FOR SLEEP
epochfilter = [];
epochfilter{1} = 'isequal($type, ''sleep'') & $sleepnum ==1';
epochfilter{2} = 'isequal($type, ''sleep'') & $sleepnum > 1';  

%Define minimum number of cells for candidate replay events
cellcountthresh = 5;

timefilter = {{'get2dstate', '((abs($velocity) < 4))'}};
clear decodefilter
global decodefilter ;
decodefilter = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
decodefilter = setfilterfunction(decodefilter, 'getpopulationevents_nolinpos', {'spikes','linpos', 'pos','ripples','cellinfo'},cellcountthresh);
decodefilter = runfilter(decodefilter);


%% DETERMINE WHICH EVENTS ARE REPLAY EVENTS
load '/data13/mcarr/RipplePaper/decodefilterB.mat'
decodefilter = decodefilterB; clear decodefilterB
for an = 1:length(decodefilter)
    for d = 1:length(decodefilter(an).output)
        for e = 1:length(decodefilter(an).output{d})
            %Determine which training index to use
            if ~isempty(decodefilter(an).output{d}(e).index)
                decode_day = decodefilter(an).output{d}(e).index(1,1);
                training_e = find(decode_day==trainingfilter(an).epochs{1}(:,1));
                
                out = calcepochreplaystats([an d e], [1 training_e], trainingfilter, decodefilter);
                if ~isempty(out)
                    ind = zeros(length(decodefilter(an).output{d}(e).eventtraj),1);
                    eventind = out(:,5);
                    decodefilter(an).output{d}(e).pvalue = NaN(size(ind));
                    decodefilter(an).output{d}(e).pvalue(eventind) = out(:,3);
                    decodefilter(an).output{d}(e).rvalue = NaN(size(ind));
                    decodefilter(an).output{d}(e).rvalue(eventind) = out(:,2);
                    decodefilter(an).output{d}(e).replaylength = NaN(size(ind));
                    decodefilter(an).output{d}(e).replaylength(eventind) = out(:,4);
                end
            end
        end
    end
end
decodefilterB = decodefilter; clear decodefilter
save('/data13/mcarr/RipplePaper/decodefilterB.mat','decodefilterB');

load '/data13/mcarr/RipplePaper/decodefilterA.mat'
decodefilter = decodefilterA; clear decodefilterA decodefilterB
for an = 1:length(decodefilter)
    for d = 1:length(decodefilter(an).output)
        for e = 1:length(decodefilter(an).output{d})
            if ~isempty(decodefilter(an).output{d}(e).index(1,1);
                decode_day = decodefilter(an).output{d}(e).index(1,1);
                training_e = find(decodeday==trainingfilter(an).epochs{2}(:,1));
                out = calcepochreplaystats([an d e], [2 training_e], trainingfilter, decodefilter);
                if ~isempty(out)
                    ind = zeros(length(decodefilter(an).output{d}(e).eventtraj),1);
                    eventind = out(:,5);
                    decodefilter(an).output{d}(e).pvalue = NaN(size(ind));
                    decodefilter(an).output{d}(e).pvalue(eventind) = out(:,3);
                    decodefilter(an).output{d}(e).rvalue = NaN(size(ind));
                    decodefilter(an).output{d}(e).rvalue(eventind) = out(:,2);
                    decodefilter(an).output{d}(e).replaylength = NaN(size(ind));
                    decodefilter(an).output{d}(e).replaylength(eventind) = out(:,4);
                end
            end
        end
    end
end
decodefilterA = decodefilter;

save('/data13/mcarr/RipplePaper/trainingfilter.mat','trainingfilter');
save('/data13/mcarr/RipplePaper/decodefilterA.mat','decodefilterA');
