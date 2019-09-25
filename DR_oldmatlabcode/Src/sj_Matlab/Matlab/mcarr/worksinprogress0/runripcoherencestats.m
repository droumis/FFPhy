%% PLOT COHERENCE

%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';

%Define iterator
iterator = 'epocheegnonreferenceanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

%Run filter for CA1 across hemispheres
CA1tetfilter = 'isequal($area,''CA1'') & $numcells> 0 & isequal($hemisphere,''right'')';
fCA1 = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',CA1tetfilter,'excludetime',timefilter,'iterator', iterator);
fCA1 = setfilterfunction(fCA1, 'loadriptriggeredcoherence', {'coherence','tetinfo'},'pair','isequal($hemisphere,''left'') & isequal($area,''CA1'') & $numcells>0');
fCA1 = runfilter(fCA1);

%% Plot average coherence
C11 = [];
for an = 1:length(fCA1)
    for d = 1:length(fCA1(an).output)
        for e = 1:length(fCA1(an).output{d})
            if isempty(C11)
                C11 = fCA1(an).output{d}(e).coherence;
            else
                C11 = cat(3,C11,fCA1(an).output{d}(e).coherence);
            end
        end
    end
end

figure
surf(time,freq,squeeze(mean(C11,3))');
view(0,90);
shading interp
set(gca,'xlim',[time(1) time(end)],'ylim',[freq(1) freq(end)])
set(gca,'clim',[0.4 0.8])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')

%Look at phase
time = fCA1(1).output{1}(1).time;
bin = -pi+pi/40:pi/20:pi-pi/40;
gam = [];
for an = 1:length(fCA1)
    animal_count = 0;
    for d = 1:length(fCA1(an).eegdata{1})
        tmp_gam = [];
        for e = 1:length(fCA1(an).eegdata{1}{d})
            animal_count = animal_count+1;
            if isempty(tmp_gam)
                tmp_gam = fCA1(an).output{1}(animal_count).gamma_phase;
                tet_count = 1;
            else
                tmp_gam = tmp_gam+fCA1(an).output{1}(animal_count).gamma_phase;
                tet_count = tet_count+1;
            end
        end
        if isempty(gam)
            gam = tmp_gam./tet_count;
        elseif ~isempty(tmp_gam)
            gam = cat(3,gam,tmp_gam./tet_count);
        end
    end
end

figure
surf(bin,time,mean(gam,3));
view(0,90);
shading interp
set(gca,'clim',[0 0.2])
set(gca,'xlim',[bin(2) bin(end-1)],'ylim',[time(1) time(end)])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')


%% Run filter for CA1-CA3 across hemispheres
CA13tetfilter = 'isequal($area,''CA1'') & $numcells> 1';
fCA13 = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',CA13tetfilter,'excludetime',timefilter,'iterator', iterator);
fCA13 = setfilterfunction(fCA13, 'loadriptriggeredcoherence', {'coherence','tetinfo'},'pair','isequal($area,''CA3'') & $numcells>1');
fCA13 = runfilter(fCA13);

%% Plot average coherence
C13 = [];
for an = 1:length(fCA13)
    for d = 1:length(fCA13(an).output)
        for e = 1:length(fCA13(an).output{d})
            if isempty(C13)
                time = fCA13(an).output{d}(e).time;
                freq = fCA13(an).output{d}(e).frequency;
                C13 = fCA13(an).output{d}(e).coherence;
            else
                C13 = cat(3,C13,fCA13(an).output{d}(e).coherence);
            end
        end
    end
end
C13 = squeeze(mean(C13,3));
figure
surf(time,freq,C13');
view(0,90);
shading interp
set(gca,'clim',[0.4 0.8])
set(gca,'xlim',[time(1) time(end)],'ylim',[freq(1) freq(end)])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')

%Look at phase
time = fCA13(1).output{1}(1).time;
bin = -pi+pi/40:pi/20:pi-pi/40;
gam = [];
for an = 1:length(fCA13)
    animal_count = 0;
    for d = 1:length(fCA13(an).eegdata{1})
        tmp_gam = [];
        for e = 1:length(fCA13(an).eegdata{1}{d})
            animal_count = animal_count+1;
            if isempty(tmp_gam)
                tmp_gam = fCA13(an).output{1}(animal_count).gamma_phase;
                tet_count = 1;
            else
                tmp_gam = tmp_gam+fCA13(an).output{1}(animal_count).gamma_phase;
                tet_count = tet_count+1;
            end
        end
        if isempty(gam)
            gam = tmp_gam./tet_count;
        else
            gam = cat(3,gam,tmp_gam./tet_count);
        end
    end
end

figure
surf(bin,time,mean(gam,3));
view(0,90);
shading interp
set(gca,'clim',[0 0.2])
set(gca,'xlim',[bin(2) bin(end-1)],'ylim',[time(1) time(end)])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')

%% Run filter for CA3-CA3 across hemispheres
CA33tetfilter = 'isequal($area,''CA3'') & $numcells> 1 & isequal($hemisphere,''right'')';
fCA33 = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',CA33tetfilter,'excludetime',timefilter,'iterator', iterator);
fCA33 = setfilterfunction(fCA33, 'loadriptriggeredcoherence', {'coherence','tetinfo'},'pair','isequal($hemisphere,''left'') & isequal($area,''CA3'') & $numcells>1');
fCA33 = runfilter(fCA33);

% Plot average coherence
C33 = [];
for an = 1:length(fCA33)
    for d = 1:length(fCA33(an).output)
        for e = 1:length(fCA33(an).output{d})
            if isempty(C33)
                time = fCA33(an).output{d}(e).time;
                freq = fCA33(an).output{d}(e).frequency;
                C33 = fCA33(an).output{d}(e).coherence;
            else
                C33 = cat(3,C33,fCA33(an).output{d}(e).coherence);
            end
        end
    end
end
C33 = squeeze(mean(C33,3));
figure
surf(time,freq,C33');
view(0,90);
shading interp
set(gca,'clim',[0.4 0.8])
set(gca,'xlim',[time(1) time(end)],'ylim',[freq(1) freq(end)])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')

%Look at phase
time = fCA33(1).output{1}(1).time;
bin = -pi+pi/40:pi/20:pi-pi/40;
gam = [];
for an = 1:length(fCA33)
    animal_count = 0;
    for d = 1:length(fCA33(an).eegdata{1})
        tmp_gam = [];
        for e = 1:length(fCA33(an).eegdata{1}{d})
            animal_count = animal_count+1;
            if isempty(tmp_gam)
                tmp_gam = fCA33(an).output{1}(animal_count).gamma_phase;
                tet_count = 1;
            else
                tmp_gam = tmp_gam+fCA33(an).output{1}(animal_count).gamma_phase;
                tet_count = tet_count+1;
            end
        end
        if isempty(gam)
            gam = tmp_gam./tet_count;
        else
            gam = cat(3,gam,tmp_gam./tet_count);
        end
    end
end

figure
surf(bin,time,mean(gam,3));
view(0,90);
shading interp
set(gca,'clim',[0 0.3])
set(gca,'xlim',[bin(2) bin(end-1)],'ylim',[time(1) time(end)])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')
