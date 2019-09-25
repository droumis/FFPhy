%% RUN FILTER AND LOOK ACROSS DAYS

%animal selection
animals = {'Conley','Corriander','Dudley','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilterA = cell(10,1);

for i = 1:10
    epochfilterA{i} = ['($exposure == ',num2str(i),') & $dailytracksexperienced == 1 & isequal($description,''TrackA'')'];
end

% Tetrode selection
ca1tetfilter =  '(isequal($area, ''CA1'') & $numcells>2 )';

% Time selection
timefilter = {};

%Select iterator
iterator = 'epocheeganal';

% Create and Run Filter
A = createfilter('animal',animals,'epochs',epochfilterA,'excludetime',timefilter, 'eegtetrodes',ca1tetfilter,'iterator', iterator);
A = setfilterfunction(A,'calceegpowerspeed',{'eeg','pos'});
A = runfilter(A);

save('/data13/mcarr/VelocityPaper/Theta/acrossdays.mat','A')

%% Look at correlation between speed and theta power and theta frequency
rho = cell(length(A),1);
powerboot = cell(length(A),1);
freqboot = cell(length(A),1);
nboot = 1000;
for an = 1:length(A)
    rho{an} = nan(length(A(an).output),4);
    powerboot{an} = nan(length(A(an).output),nboot);
    freqboot{an} = nan(length(A(an).output),nboot);
    for d = 1:length(A(an).output)
        if ~isempty(A(an).output{d})
        velocity = A(an).output{d}(1).speed;
        tmp = nan(length(A(an).output{d}(1).power),length(A(an).output{d}));
        for e = 1:length(A(an).output{d})
            tmp(:,e) = A(an).output{d}(e).power;
        end
        power = nanmean(tmp,2);
        tmp = nan(length(A(an).output{d}(1).frequency),length(A(an).output{d}));
        for e = 1:length(A(an).output{d})
            tmp(:,e) = A(an).output{d}(e).frequency;
        end
        frequency = nanmean(tmp,2);
        power(isnan(velocity) | isinf(velocity)) = [];
        frequency(isnan(velocity) | isinf(velocity)) = [];
        velocity(isnan(velocity) | isinf(velocity)) = [];
        [r1 p1] = corr(velocity,power);
        [r2 p2] = corr(velocity,frequency);
        rho{an}(d,:) = [r1 p1 r2 p2];
        end
        
        for s = 1:nboot
            boot = ceil(length(velocity)*rand(length(velocity),1));
            xboot = velocity(boot);
            y1boot = power(boot);
            y2boot = frequency(boot);
            b = regress(y1boot,[ones(length(xboot),1) xboot],0.05);
            powerboot{an}(d,s) = b(2);
            b = regress(y2boot,[ones(length(xboot),1) xboot],0.05);
            freqboot{an}(d,s) = b(2);
        end
    end
end

%% RUN FILTER AND LOOK AT WITHIN DAY COMPARISONS

%animal selection
animals = {'Corriander','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilterA = [];
epochfilterB = [];
for i = 1:12
    epochfilterA{i} = ['($experimentday == ',num2str(i),') & $dailyexposure == 1 & isequal($description,''TrackA'')'];
    epochfilterB{i} = ['($experimentday == ',num2str(i),') & $dailyexposure == 1 & isequal($description,''TrackB'')'];
end

% Tetrode selection
ca1tetfilter =  '(isequal($area, ''CA1'') & $numcells>2 )';

% Time selection
timefilter = {};

%Select iterator
iterator = 'epocheeganal';

% Create and Run Filter
A = createfilter('animal',animals,'epochs',epochfilterA,'excludetime',timefilter, 'eegtetrodes',ca1tetfilter,'iterator', iterator);
A = setfilterfunction(A,'calceegpowerspeed',{'eeg','pos'});
A = runfilter(A);

% Create and Run Filter
B = createfilter('animal',animals,'epochs',epochfilterB,'excludetime',timefilter, 'eegtetrodes',ca1tetfilter,'iterator', iterator);
B = setfilterfunction(B,'calceegpowerspeed',{'eeg','pos'});
B = runfilter(B);

save('/data13/mcarr/VelocityPaper/Theta/withindayA.mat','A')
save('/data13/mcarr/VelocityPaper/Theta/withindayB.mat','B')

%% Look at within day correlations for speed and theta power and frequency
rho = nan(length(A),4);
phi = nan(length(A),4);
nboot = 10000;
powerAboot = nan(length(A),nboot); powerBboot = nan(length(A),nboot);
freqAboot = nan(length(A),nboot); freqBboot = nan(length(A),nboot);

for an = 1:length(A)
    count = 0;
    for d = 1:length(A(an).output)
        if ~isempty(A(an).output{d}) && ~isempty(B(an).output{d})
            count = count+1;
            if count == 1
                
            velA = A(an).output{d}(1).speed;
            tmp = nan(length(A(an).output{d}(1).power),length(A(an).output{d}));
            tmp2 = nan(length(A(an).output{d}(1).frequency),length(A(an).output{d}));
            for e = 1:length(A(an).output{d})
                tmp(:,e) = A(an).output{d}(e).power;
                tmp2(:,e) = A(an).output{d}(e).frequency;
            end
            powA = nanmean(tmp,2);
            freqA = nanmean(tmp2,2);
            powA(isnan(velA) | isinf(velA)) = [];
            freqA(isnan(velA) | isinf(velA)) = [];
            velA(isnan(velA) | isinf(velA)) = [];
            [r p] = corr(velA,powA);
            rho(an,1:2) = [r p];
            [r p] = corr(velA,freqA);
            phi(an,1:2) = [r p];
            for s = 1:nboot
                boot = ceil(length(velA)*rand(length(velA),1));
                xboot = velA(boot);
                y1boot = powA(boot);
                y2boot = freqA(boot);
                b = regress(y1boot,[ones(length(xboot),1) xboot],0.05);
                powerAboot(an,s) = b(2);
                b = regress(y2boot,[ones(length(xboot),1) xboot],0.05);
                freqAboot(an,s) = b(2);
            end
            

            velB = B(an).output{d}(1).speed;
            tmp = nan(length(B(an).output{d}(1).power),length(B(an).output{d}));
            tmp2 = nan(length(B(an).output{d}(1).frequency),length(B(an).output{d}));
            for e = 1:length(B(an).output{d})
                tmp(:,e) = B(an).output{d}(e).power;
                tmp2(:,e) = B(an).output{d}(e).frequency;
            end
            powB = nanmean(tmp,2);
            freqB = nanmean(tmp2,2);
            powB(isnan(velB) | isinf(velB)) = [];
            freqB(isnan(velB) | isinf(velB)) = [];
            velB(isnan(velB) | isinf(velB)) = [];
            [r p] = corr(velB,powB);
            rho(an,3:4) = [r p];
            [r p] = corr(velB,freqB);
            phi(an,3:4) = [r p];
            for s = 1:nboot
                boot = ceil(length(velB)*rand(length(velB),1));
                xboot = velB(boot);
                y1boot = powB(boot);
                y2boot = freqB(boot);
                b = regress(y1boot,[ones(length(xboot),1) xboot],0.05);
                powerBboot(an,s) = b(2);
                b = regress(y2boot,[ones(length(xboot),1) xboot],0.05);
                freqBboot(an,s) = b(2);
            end
            
            end
        end
    end
end
