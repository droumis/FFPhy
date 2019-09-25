% Animal Selection
animals = {'Conley'};
%animals = {'Bond','Conley','Corriander','Dudley','Frank','Miles','Ten'};

% Epoch selection
epochfilter = [];
for i = 1
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description,''TrackA'')'];
end

% Tetrode selection
%ca1tetfilter = '(isequal($area, ''MEC''))';
tetrodefilter = '(isequal($area, ''CA1'') & ($numcells > 2))';

% Time selection
timefilter = {{'getlinstate', '($velocity >15)', 6}};

%Select iterator
iterator = 'epocheeganal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter, 'eegtetrodes',tetrodefilter,'iterator', iterator);
f = setfilterfunction(f,'calcaveragetrace',{'eeg','ripple'});
f = runfilter(f);

%% PLOT
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e))
                figure
                plot(-250/1.5:1/1.5:250/1.5,f(an).output{d}(e).eeg,'b',...
                    -250/1.5:1/1.5:250/1.5,f(an).output{d}(e).rip,'r')
            end
        end
    end
end

