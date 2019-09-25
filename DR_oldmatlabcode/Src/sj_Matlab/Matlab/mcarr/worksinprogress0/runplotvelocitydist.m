% Animal Selection
animals = {'Eight'};

% Epoch selection
epochfilter = [];
for i = 1:10
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description,''TrackA'')'];
end

% Select iterator
iterator = 'singleepochanal';

% Time selection for Low Gamma
timefilter = {{'getriptimes', '($nripples >= 1)', [], 'tetfilter', 'isequal($area,''CA1'')'}};

% Create and Run filter velocity
g = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter, 'iterator',iterator);
g = setfilterfunction(g,'plotvelocitydist',{'pos'});
g = runfilter(g);
