%% RUN FILTER FOR ALL CELLS
%Animal selection
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
%epochfilter{1} = 'isequal($type, ''sleep'')';

cellfilter = '($meanrate<7) && isequal($area,''CA1'')';


%Define iterator
iterator = 'singlecellanal';

%create training data by calulating the linearized rates of all cells
f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'iterator',iterator);
f = setfilterfunction(f, 'calcclusterquality', {'clustqual'});
f = runfilter(f);

out = numericgroupcombine(f);

% Isolation distance

% CA1: mean = 54.97, s.e.m. = 4.6, std = 165.9
% CA3: mean = 94.612, s.e.m. = 16.4190, std = 648

%L-Ratio

% CA1: mean = 23.8558, s.e.m. = 5.0764, std = 182.8
% CA3: mean = 8.73, s.e.m. = 1.33, std = 52.6
