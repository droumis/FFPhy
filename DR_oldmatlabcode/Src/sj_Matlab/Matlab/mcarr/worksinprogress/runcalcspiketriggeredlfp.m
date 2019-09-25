%% RUN FOR CA1 RIPPLES

% Animal Selection% Animal Selection
animals = {'Ten','Bond','Frank'};

% Epoch selection
epochfilter = [];
epochfilter{1} = '(isequal($type,''run''))';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

%tetrode filter
tetfilter = 'isequal($area, ''CA1'') & $numcells > 2';

iterator = 'epocheegnonreferenceanal';

f = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',tetfilter,'excludetimefilter', timefilter, 'iterator', iterator);
f = setfilterfunction(f, 'calcspiketriggeredlfp', {'ripple', 'ripples', 'spikes', 'cellinfo'});
f = runfilter(f);

save('/data21/mcarr/RipplePaper/ca1ripplespiking.mat','f')

%% RUN FOR CA3 RIPPLES

% tetrode filter
tetfilter = 'isequal($area, ''CA3'') & $numcells > 1';

iterator = 'epocheegnonreferenceanal';

f = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',tetfilter,'excludetimefilter', timefilter, 'iterator', iterator);
f = setfilterfunction(f, 'calcspiketriggeredlfp', {'ripple', 'ripples', 'spikes', 'cellinfo'});
f = runfilter(f);

save('/data21/mcarr/RipplePaper/ca3ripplespiking.mat','f')

%% PLOT

x = []; s = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            for c= 1:length(f(an).output{d}(e).lfp)
                x = [x; f(an).output{d}(e).lfp{c}];
                s = [s f(an).output{d}(e).ripstd{c}];
            end
        end
    end
end

%Try estimating the phase at each spike
time = -0.01:1/1500:0.01;
zero_time = lookup(0,time);
for c = 1:size(x,1)
    tmp = x(c,:);
    
    %Identify the zero crossing before and after
    zero_before = [find(tmp(1:zero_time)>0,1,'last') find(tmp(1:zero_time)<0,1,'last')];
    zero_after = [zero_time+find(tmp(zero_time:end)<0,1,'first') zero_time+find(tmp(zero_time:end)>0,1,'first')];
    
    [maxtab mintab] = peakdet(tmp,0.5);
    previous_peak = find(maxtab(:,1)<zero_time,1,'last');
    next_peak = find(maxtab(:,1)>= zero_time,1,'first');
    
    %Create linear spaced phase between these two points
    tmp_time = linspace(time(maxtab(previous_peak,1)),time(maxtab(next_peak,1)),360);
    phase(c) = lookup(0,tmp_time);
end
tmp = phase;
tmp(tmp>180) = tmp(tmp>180)-360;

%Plot spike triggered ripple trace for each threshold
threshbin =[0 3.5; 0 5; 0 6; 0 7; 0 100];

R1 = zeros(length(time),size(threshbin,1));
P1 = zeros(17,size(threshbin,1));
count = 1;
for i = threshbin'
    R1(:,count) = median(x(s>i(1) & s<i(2),:));
    P1(:,count) = hist([tmp(s>i(1) & s<i(2)) 360+tmp(s>i(1) & s<i(2))],-180:45:540);
    count = count+1;
end

color = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1];
figure
for i = 1:size(R1,2)
    plot(time,R1(:,i),'color',color(i,:))
    hold on
end
legend([{'CA1 < 4std'},{'CA1<5std'},{'CA1<6std'},{'CA1<7std'},{'CA1<20std'}])
xlabel('Time relative to spike (s)')
ylabel('Average CA3 Ripple Filtered LFP (150-250HZ)')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_CA3_spiketriggeredripplelfp.pdf', m, d, y);
print('-dpdf', savestring)

