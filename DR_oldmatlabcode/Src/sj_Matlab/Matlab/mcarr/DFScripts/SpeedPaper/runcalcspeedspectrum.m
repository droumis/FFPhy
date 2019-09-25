%% RUN FILTER FOR CA1 ACROSS DAYS

%animal selection
animals = {'Conley','Corriander','Dudley','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilter = [];
for i = 1
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description,''TrackA'')'];
end

% Tetrode selection
ca1tetfilter =  '(isequal($area, ''CA1'') & $numcells>=1 )';

% Time selection
timefilter = {};

%Select iterator
iterator = 'epocheeganal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter, 'eegtetrodes',ca1tetfilter,'iterator', iterator);
f = setfilterfunction(f,'calcspeedspectrum',{'eeg','pos'});
f = runfilter(f);

speed_spectrum = f;
save('/data13/mcarr/VelocityPaper/speedspectrum.mat','speed_spectrum')

%% LOad
load('/data13/mcarr/VelocityPaper/speedspectrum.mat')
f = speed_spectrum;
clear speed_spectrum

%% PLOT
for an = 2%1:length(f)
    for day = 1%1:length(f(an).output)
        freq = f(an).output{day}(1).frequency;
        bin = f(an).output{day}(1).bin;
        tmp = zeros(size(f(an).output{day}(1).spectrum));
        for tet = 1:length(f(an).output{day})
            tmp = [tmp + f(an).output{day}(tet).spectrum];
        end
        tmp = tmp./length(f(an).output{day});
        figure
        surface(mean(bin,2), freq, tmp,'EdgeColor','none','FaceColor','interp')
        set(gca,'clim',[-0.75 0.75])
        colorbar
        set(gca,'xscale','log')
        set(gca,'xlim',[0 max(mean(bin,2))])
        set(gca,'ylim',[freq(1) freq(end)])
        set(gca,'xtick',[1/4 1/2 1 2 4 8 16])
        set(gca,'FontSize',18)
        colormap('gray')
        %box off
    end
end

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_velocity_speed_spectrum_Fivday1.png', m, d, y);
print('-dpng', savestring)
