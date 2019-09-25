runscript = 1

%% currently don't know why not all tetrodes for chapati show up



if runscript
    

% Animal
animals = {'Chapati'};

dayfilter = 1:11;

epochfilter{1} = ['isequal($type,''run'')'];

tetfilter = '(isequal($area, ''CA3'') & $numcells >= 2)';

iterator = 'eeganal';

timefilter = { {'kk_get2dstate', '$velocity > 4'} };

f = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter','excludetime',timefilter,'iterator', iterator,'eegtetrodes',tetfilter);
f = setfilterfunction(f,'dfakk_MI',{'eeg'});
f = runfilter(f);

cfc = f.output;

end

if 0
   save('MI') 
end

%% manually rename file here


%% Plot results by tetrode across all days and epochs.

% iterate across tetrodes

% first collect all epoch indices [d e t]
indices = [];
for d = 1:length(cfc)
    for e = 1:length(cfc{d})
        indices = [indices ; cfc{d}(e).index];
    end
end
tetrodes = unique(indices(:,3))';
    
% now for each tetrode plot all epochs
for tet = tetrodes(4)
    flag = 0;
    figure
    % collect the tetrode's epochs
    tetepochs = indices(indices(:,3)==tet,:);
    % sort epochs chronologically
    tetepochs = sortrows(tetepochs,[1 2]);
    numepochs = size(tetepochs,1);
    for ep = 1:numepochs
        figure
        %subplot(1,numepochs,ep)
        ind = rowfind(tetepochs(ep,:),indices);
%        imagesc(cfc{d}(ind).frequency_phase, ...
%                cfc{d}(ind).frequency_amplitude, ...
%                cfc{d}(ind).comodulogram')

            contourf(cfc{d}(ind).frequency_phase +cfc{d}(end).phasefreq_bandwidth/2,...
                 cfc{d}(ind).frequency_amplitude+cfc{d}(end).ampfreq_bandwidth/2, ...
                 cfc{d}(ind).comodulogram',30,'lines','none')
              colorbar
            %caxis([.00001 .001]);

        set(gca,'YDir','normal');
        colormap('hot')
        title([num2str(tetepochs(ep,1)) '  ' num2str(tetepochs(ep,2)) ' ' num2str(tetepochs(ep,3))])
    end
end





