%% Chapati
for i = 1:1
% Animal selection
animals = {'Chapati'};

% Epoch selection
epochfilter = [];
epochfilter{1} = ['isequal($type,''run'')'];

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
cfc = loaddatastruct('/datatmp/kkay/ProcessedData/','cha','cfc');
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