% Animal selection
animals = {'Five'};

% Epoch selection
epochfilter = [];
%epochfilter{1} = ['(isequal($type, ''run''))'];
epochfilter{1} = 'isequal($description,''TrackB'') & ($exposure == 1)';

% Tetrode selection
tetrodefilter = '(isequal($area, ''CA1'') & $numcells>2)';

% Time selection
timefilter = {{'getlinstate', '(abs($velocity) >3)', 6}};

% Iterator selection
iterator = 'eeganal';


% Create and Run Filter
%f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodes', tetrodefilter, 'iterator', iterator);
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodes', tetrodefilter, 'iterator', iterator);

%f = setfilterfunction(f, 'plotspectrogram', {'eeg'},'appendindex',1,'fpass',[1 200]);
f = setfilterfunction(f, 'calcspectrogram', {'eeg'},'appendindex',1,'fpass',[1 300],'window',[0.5 0.5]);

f = runfilter(f);

%% Plot the Spectrogram

for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            figure
            imagesc(f(an).output{d}(e).time,f(an).output{d}(e).frequency, ...
                log10(f(an).output{d}(e).fullspectrum'))
            axis xy
            hold on
        end
    end
end

%% Try looking for defining characteristics of ripples

f1 = [150 250];
f2 = [70 250; 150 300];

for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            freq1 = lookup(f1,f(an).output{d}(e).frequency);
            power1 = nan(size(f(an).output{d}(e).fullspectrum(:,freq1(1):freq1(2))));
            for s = freq1(1):freq1(2)
                power1(:,s+1-freq1(1)) = f(an).output{d}(e).fullspectrum(:,s)./median(f(an).output{d}(e).fullspectrum(:,s));
            end
            power1 = sum(power1,2);
            freq2 = lookup(f2,f(an).output{d}(e).frequency);
            power2 = nan(size(f(an).output{d}(e).fullspectrum(:,[freq2(1):freq2(2) freq2(3):freq2(4)])));
            for s = [freq2(1):freq2(2) freq2(3):freq2(4)]
                power2(:,s+1-freq2(1)) = f(an).output{d}(e).fullspectrum(:,s)./median(f(an).output{d}(e).fullspectrum(:,s));
            end
            power2 = nansum(power2,2);
            figure
            imagesc(f(an).output{d}(e).time,f(an).output{d}(e).frequency, ...
                log10(f(an).output{d}(e).fullspectrum'))
            axis xy
            hold on
            plot(f(an).output{d}(e).time,log10(power1)./log10(power2));
        end
    end
end

%% Plot the Spectrogram

for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            spec = nan(size(f(an).output{d}(e).fullspectrum));
            for s = 1:length(f(an).output{d}(e).frequency)
                spec(:,s) = f(an).output{d}(e).fullspectrum(:,s)./...
                    median(f(an).output{d}(e).fullspectrum(:,s));
            end
            figure
            imagesc(f(an).output{d}(e).time,f(an).output{d}(e).frequency,spec')
            axis xy
            hold on
        end
    end
end
