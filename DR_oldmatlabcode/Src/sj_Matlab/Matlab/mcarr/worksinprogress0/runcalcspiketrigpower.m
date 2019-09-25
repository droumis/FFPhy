% Animal selection
animals = {'Five'};

% Epoch selection
epochfilter = [];
epochfilter{1} = ['isequal($description, ''TrackA'') & ($exposure == 1)'];

% Tetrode selection
ca1tetfilter = '(isequal($area, ''CA1'') & ($representative ==1))';

% Cell Selection
ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 7) && ($numspikes > 200))';

% Time selection
timefilter = {{'getlinstate', '(($traj ~= -1) & abs($velocity) > 10)', 6}};

% Iterator selection
iterator = 'singlecelleeganal';


% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'cells', ca1cellfilter, 'eegtetrodes',ca1tetfilter,'iterator', iterator);

f = setfilterfunction(f, 'calcspiketrigpower', {'eeg','spikes'},'appendindex',1,'fpass',[1 250]);

f = runfilter(f);

%% TRY TO COMBINE CELLS FOR A GIVEN DAY

for an = 1:length(f)
    for d = 1:length(f(an).epochs)
        for e = 1:length(f(an).output{d})
            temp = mean(f(an).output{d}(e).power,3);
            if isempty(avg)
                avg = temp;
            else
                avg = avg + temp;
            end
            imagesc(f(an).output{d}(e).time,f(an).output{d}(e).frequency,temp')
            axis xy
            set(gca,'clim',[0 10]);
            pause
        end
    end
end