%% BOND
for i = 7:14
% Animal selection
animals = {'Bond'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs

for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;       
    end
save(cfcfile, 'cfc');
end
end

%% CONLEY
for i = 1:6
% Animal selection
animals = {'Conley'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

% Save filter outputs

for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;       
    end
save(cfcfile, 'cfc');
end
end

%% CORRIANDER
for i = 1:9
% Animal selection
animals = {'Corriander'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs

for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;       
    end
save(cfcfile, 'cfc');
end
end

%% DUDLEY
for i = 1:6
% Animal selection
animals = {'Dudley'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs
for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;       
    end
save(cfcfile, 'cfc');
end
 
end

%% Eight
for i = 1:7
% Animal selection
animals = {'Eight'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type,''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs
for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;
    end
save(cfcfile, 'cfc');
end
 
end

%% Five
for i = 1:9
% Animal selection
animals = {'Five'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs

for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;       
    end
save(cfcfile, 'cfc');
end
end

%% Frank
for i = 7:17
% Animal selection
animals = {'Frank'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs

for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;       
    end
save(cfcfile, 'cfc');
end
end

%% Miles
for i = 1:5
% Animal selection
animals = {'Miles'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs
for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;
    end
save(cfcfile, 'cfc');
end
 
end

%% SEVEN
for i = 1:9
% Animal selection
animals = {'Seven'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs

for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;       
    end
save(cfcfile, 'cfc');
end
end

%% SIX
for i = 1:5
% Animal selection
animals = {'Six'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs

for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;       
    end
save(cfcfile, 'cfc');
end
end


%% TEN
for i = 1:7
% Animal selection
animals = {'Ten'};

% Epoch selection% Epoch selection
epochfilter = [];
epochfilter{1} = ['($experimentday == ',num2str(i),') & isequal($type, ''run'')'];

% Tetrode selection
tetfilter = '(isequal($area, ''CA1'') & $numcells>2)';

%Iterator Selection
iterator = 'eeganal';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

f = createfilter('animal',animals,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);

%% Save filter outputs
for d=1:length(f.output)
    cfc = {};
    for e=1:length(f.output{d})
        if f.output{d}(e).index(2)<10
            day = sprintf('0%s',num2str(f.output{d}(e).index(1)));
        else
            day = num2str(f.output{d}(e).index(1));
        end
        epoch = f.output{d}(e).index(2);
        tet = f.output{d}(e).index(3);
        cfcfile = sprintf('%s%scfc%s.mat', ...
            f.animal{2}, f.animal{3}, day);
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.cfc = f.output{d}(e).cfc;
        cfc{f.output{d}(e).index(1)}{epoch}{tet}.frequency = f.output{d}(e).frequency;       
    end
save(cfcfile, 'cfc');
end
 
end

%% Plot results
cfc = loaddatastruct('/data13/mcarr/Eig/','Eig','cfc');
for day = 1:length(cfc);
    tmp = [];
    for epoch = 2:2:length(cfc{day})
        for tet = 1:length(cfc{day}{epoch})
            if ~isempty(cfc{day}{epoch}{tet})
                if isempty(tmp)
                    tmp = cfc{day}{epoch}{tet}.cfc - ...
                        min(min(cfc{day}{epoch}{tet}.cfc));
                    freq = cfc{day}{epoch}{tet}.frequency;
                    count = 1;
                else
                    tmp = tmp + cfc{day}{epoch}{tet}.cfc ...
                        - min(min(cfc{day}{epoch}{tet}.cfc));
                    count = count+1;
                end
            end
        end
    end
    tmp = tmp./count;
    figure
    subplot(1,2,1)
    imagesc(freq,20:2:200,tmp)
    axis xy
    subplot(1,2,2)
    plot(20:2:200,mean(tmp(:,7:13),2))
end