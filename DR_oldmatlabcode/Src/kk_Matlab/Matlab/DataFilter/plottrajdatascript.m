Veqn = '>=4'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 3 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6


%Animal selection
%-----------------------------------------------------
animals = {'Barack', 'Calvin', 'Dwight'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------



%Filter creation
%--------------------------------------------------------
epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
epochfPF = ['($switchday > 0)'];

cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 4))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal

timefilter = { {'getlinvelocity', ['((abs($velocity)', Veqn,'))']} };

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);

%only include cells with placefields
if minPeakPF>0
    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
    f = excludecellsf(f, includecells);
end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';

f = setfilteriterator(f,iterator);

f = setfilterfunction(f, 'plottrajdata', {'spikes', 'linpos'}, includestates, minV);
% out = plottrajdata(index, excludetimes, spikes, linpos, includestates, minV, varargin)

f = runfilter(f);
