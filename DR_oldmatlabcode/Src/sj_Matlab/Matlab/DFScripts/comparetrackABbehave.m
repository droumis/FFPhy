%Animal selection
%-----------------------------------------------------
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
animals = {'Conley','Bond','Frank','Ten','Alex'};
%animals = {'Ten'};
%aniamls = {'Alex'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Bond'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];

epochfilter{1} = ['(($exposureday == 1) & ($dailyexposure == 1) & isequal($environment, ''TrackB''))'];
epochfilter{2} = ['(($exposureday == 1) & ($dailyexposure == 2) & isequal($environment, ''TrackB''))'];
epochfilter{3} = ['(($exposureday == 2) & ($dailyexposure == 1) & isequal($environment, ''TrackB''))'];
epochfilter{4} = ['(($exposureday == 2) & ($dailyexposure == 2) & isequal($environment, ''TrackB''))'];
epochfilter{5} = ['(($exposureday == 3) & ($dailyexposure == 1) & isequal($environment, ''TrackB''))'];
epochfilter{6} = ['(($exposureday == 3) & ($dailyexposure == 2) & isequal($environment, ''TrackB''))'];

%cellfilter = [];


%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 3))', 6} };
timefilter = {};

%f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter);
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime', timefilter);


%-----------------------------------------------------------




%run function- single cells
%--------------------------------------------
iterator = 'epochbehaveanal';

f = setfilteriterator(f,iterator);
f = setfilterfunction(f, 'calcPassTimes', {'linpos'});
f = runfilter(f);

for a = 1:5;
    
    day1 = [f(a).output{1}.passtimes;f(a).output{2}.passtimes];
    day2 = [f(a).output{3}.passtimes;f(a).output{4}.passtimes];
    day3 = [f(a).output{5}.passtimes;f(a).output{6}.passtimes];
    day1mean(a) = mean(day1);
    day2mean(a) = mean(day2);
    day3mean(a) = mean(day3);
    day1stderr(a) = stderr(day1);
    day2stderr(a) = stderr(day2);
    day3stderr(a) = stderr(day3);
    
end

figure
bar(1:5,day1mean(1:5))
hold on
errorbar(1:5,day1mean(1:5),day1stderr(1:5),'.')
bar([1:5],day2mean(1:5),'r')
errorbar(1:5,day2mean(1:5),day2stderr(1:5),'.')

bar([1:5],day3mean(1:5),'g')
errorbar(1:5,day3mean(1:5),day3stderr(1:5),'.')

    



