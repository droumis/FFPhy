runscript = 1

%% currently don't know why not all tetrodes for chapati show up



if runscript
    

% Animal
animals = {'Egypt'};

dayfilter = 1:11;

epochfilter{1} = ['isequal($type,''run'')'];

tetfilter = '(isequal($area, ''CA3'') & $numcells > 2)';

iterator = 'eeganal';

timefilter = { {'kk_get2dstate', '$velocity > 4'} };

f = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'crossfrequencycoherence',{'eeg'});
f = runfilter(f);


if 1
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
end

%% Plot results by tetrode across all days and epochs!

% iterate across tetrodes

% first collect all epochs each tetrode participates in
plotmatrix = {};
for tet = 1:21
    plotmatrix{tet} = [];
    for d = 1:length(cfc)
        for e = 1:length(cfc{d})
            if ~isempty(cfc{d}{e})
                try
                    if ~isempty(cfc{d}{e}{tet})
                        plotmatrix{tet} = [plotmatrix{tet} ; d e];
                    end
                catch
                end
            end
        end
    end 
end

% now for each tetrode plot all epochs
for t = 1:length(plotmatrix)
    if ~isempty(plotmatrix{t})
        figure
        for k = 1:size(plotmatrix{t},1)
            subplot(1,size(plotmatrix{t},1),k)
            day = plotmatrix{t}(k,1);
            epoch = plotmatrix{t}(k,2);
            imagesc(cfc{day}{epoch}{t}.frequency(5:18),20:2:200,cfc{day}{epoch}{t}.cfc(:,5:18),[.425 .47])
                set(gca,'YDir','normal');
                colormap('hot')      
            title([num2str(t) '  ' num2str(day) ' ' num2str(epoch)])
        end
    end
end




