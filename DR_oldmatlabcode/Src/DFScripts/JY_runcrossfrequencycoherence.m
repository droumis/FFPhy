%% BOND

days='[5]';%,'1:10';

% Animal selection
animals = {'N2'};

% Epoch selection% Epoch selection
%epochfilter = [];

%epochfilter = ['isequal($epochtype, ''Run'')'];
epochfilter = ['isequal($epoch,  6)'];

% Tetrode selection
tetfilter = '(isequal($area, ''ACC''))';

%Iterator Selection
iterator = 'JY_eeganal';

% Time selection
timefilter = {{'JY_getlinvelocity', 'abs($velocity) < 3'}};
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'JY_crossfrequencycoherence',{'eeg','meaneegspectrograms'});
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
%save(cfcfile, 'cfc');
end

% quick plot

cfcout=arrayfun(@(x) x.cfc,f.output{1,1},'Uniformoutput',false);
freqsx=arrayfun(@(x) x.frequency,f.output{1,1},'Uniformoutput',false);
freqsy=arrayfun(@(x) x.amp_freq,f.output{1,1},'Uniformoutput',false);
for ii=1:size(cfcout,2)
    tmp=cfcout{ii};
    freqx=freqsx{ii};
    freqy=freqsy{ii};
    figure
    %subplot(1,2,1)
    imagesc(freqx,freqy,tmp)
    colorbar
    axis xy
    %subplot(1,2,2)
    %plot(20:2:200,mean(tmp,1))
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