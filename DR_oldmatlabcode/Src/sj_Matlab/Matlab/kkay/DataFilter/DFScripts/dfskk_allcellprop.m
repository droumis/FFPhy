
%Animal selection
%-----------------------------------------------------
%animals = {'Bond','Frank','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Frank'};
animals = {'Egypt'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = 'isequal($type, ''run'')';
%epochfilter{1} = ['isequal($environment, ''TrackB'') && ($dailyexposure == 1)'];
%for i = 1:14
    %epochfilter{i} = ['isequal($environment, ''TrackA'') | isequal($environment, ''TrackB'')'];
%   epochfilter{i} = ['isequal($environment, ''TrackA'')'];
%nd

ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100))';

% 1. nothing
%timefilter = {};

% 2. moving around, no ripples
timefilter = { {'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...     % exclude ripples
               {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 5))', 6} };                            % velocity cutoff

% 3. ripple
%timefilter = {{'getriptimes', '($nripples > 1)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};


iterator = 'singlecellanal';
%iterator = 'multicellanal';

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'iterator', iterator, 'excludetime', timefilter);
ca2f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca2cellfilter,'iterator', iterator, 'excludetime', timefilter);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'iterator', iterator, 'excludetime', timefilter);
%-----------------------------------------------------------


ca1f = setfilterfunction(ca1f, 'cellprop', {'spikes'});
ca2f = setfilterfunction(ca2f, 'cellprop', {'spikes'});
ca3f = setfilterfunction(ca3f, 'cellprop', {'spikes'});


'ca1'
ca1f = runfilter(ca1f);
ca1 = numericgroupcombine(ca1f);

'ca2'
ca2f = runfilter(ca2f);
ca2 = numericgroupcombine(ca2f);

'ca3'
ca3f = runfilter(ca3f);
ca3 = numericgroupcombine(ca3f);


% mean_waveform_width  mean_rate  csi  prop_bursts *theta_mod_depth  
%   *peak_theta_phase gamma_mod_depth  peak_gamma_phase

for i=1:4
figure; hold on;
[N X] = hist(ca1{1}(:,i));
hist(ca1{1}(:,i));
hist(ca2{1}(:,i),X);
hist(ca3{1}(:,i),X)
h = findobj(gca,'Type','patch');
display(h)
set(h(1),'FaceColor','g','EdgeColor','w')
alpha(0.2)
set(h(2),'FaceColor','k','EdgeColor','w')
alpha(0.2)
set(h(3),'FaceColor','r','EdgeColor','w')
alpha(0.2)
end















