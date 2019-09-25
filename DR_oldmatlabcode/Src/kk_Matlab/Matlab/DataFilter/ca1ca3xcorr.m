
%Animal selection
%-----------------------------------------------------
animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];

cellpairfilter = {'allcomb', '(isequal($area, ''CA3'') && ($meanrate < 4))', '(isequal($area, ''CA3'') && ($meanrate < 4))'};
%makes pairs, always takes 3 arguments, first argument is how to combine
%other 2 arguments, arg 2 is 1st member of pair, 2nd is 2nd member of pair

timefilter = {{'getriptimes', '($nripples >= 2)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};
%nripp[les >=2 means ripple on at least 2 tetrodes, need to give it a lits
%of tetrodes for that: see help getriptimes
%plot distribution of energies 

iterator = 'singlecellanal';
%each pair at once

f = createfilter('animal',animals,'epochs',epochfilter,'cellpairs',cellpairfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'calcxcorrmeasures', {'spikes'}, 'edgespikes', 1);
f = runfilter(f);

%calcxcorrmeasures does cross correlation stuff, taskes cell pairs in
%to do ripple rate analysis use cellfilter from prev, need new function &
%new iterator
