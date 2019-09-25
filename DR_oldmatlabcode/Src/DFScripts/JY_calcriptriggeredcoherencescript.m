%% Run Filter
animals = {'N2'};
%MECanimals = {'Eight','Corriander'};
days='[5]';
% Epoch selection
epochfilter = [];
% for i = 1:24
%     epochfilter{i} = ['($exposure == ',num2str(i),')'];
% end
epochfilter = ['isequal($epoch,  4)'];

% Tetrode selection
CA3tetrodepairfilter = {'(isequal($area, ''ACC''))', '(isequal($area, ''HP'')) '};
%MECtetrodepairfilter = {'(isequal($area, ''CA1'') & ($numcells > 2))', '(isequal($area, ''MEC''))'};

% Time selection
timefilter = {{'JY_getlinvelocity', 'abs($velocity) >3'}};

% Iterator selection
iterator = 'JY_eegnonreferenceanal';

% Create and Run Filter
%f = createfilter('animal',MECanimals,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodepairs', MECtetrodepairfilter, 'iterator', iterator);
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'excludetime',timefilter,'eegtetrodepairs', CA3tetrodepairfilter, 'iterator', iterator);

f = setfilterfunction(f, 'JY_calcriptriggeredcoherence', {'data','meaneegspectrograms','ripples'},'rippletetrode',[18],'appendindex',1);

f = runfilter(f);

%% Load Data

% concatenate spectra

Cconc=arrayfun(@(x) x.coherence,f.output{1,1},'uniformoutput',0);
Cconc=cat(3,Cconc{:});
Cconcmean=mean(Cconc,3);

t=f.output{1,1}(1,1).times;
freq=f.output{1,1}(1,1).frequency;
halfwindow=2;
figure;

imagesc(t-halfwindow/10000,freq,Cconcmean');
axis xy;
colorbar;

% Load data and define options
load '/data13/mcarr/VelocityPaper/Coherence/ca3ca1nonreferencecoherence.mat'
c = ca3ca1nonreferencecoherence;
clear ca3ca1nonreferencecoherence

load '/data13/mcarr/VelocityPaper/Coherence/mecca1nonreferencecoherence.mat'
m = mecca1nonreferencecoherence;
clear mecca1nonreferencecoherence

frequency = c(1).output{end}(end).frequency;
noise = [lookup(59,frequency) lookup(61,frequency)];
range = lookup(140,frequency);
freq = frequency([1:noise(1)-1 noise(2)+1:range]);
g = gaussian(5,25);
%% Plot examples

% Plot example Coherence for CA3-CA1
call = [];
for an = 3
    for d = 1:length(c(an).output)
        tmp = [];
        for e = 1:length(c(an).output{d})
            if ~isempty(c(an).output{d})
                 tmp = [tmp; ...
                    c(an).output{d}(e).coherence([1:noise(1)-1 noise(2)+1:range])];
            end
        end
        if size(tmp,1) > 1
             call = [call; mean(tmp)];
        end         
    end
end

figure
plot(freq,smoothvect(call(10,:),g),'k')
set(gca,'xlim',[2 139],'ylim',[0.4 0.8],'FontSize',18)
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Coherence','FontSize',18)
title('CA3 - CA1','FontSize',18)
box off

%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_ca3ca1_coherence_example.pdf', m, d, y);
print('-dpdf', savestring)

% Plot example coherence for MEC-CA1
call = [];
for an = 2
    for d = 1:length(m(an).output)
        tmp = [];
        for e = 1:length(m(an).output{d})
            if ~isempty(m(an).output{d})
                 if any(m(an).output{d}(e).index(4) == [6 7])
                    tmp = [tmp; ...
                    m(an).output{d}(e).coherence([1:noise(1)-1 noise(2)+1:range])];
                 end
            end
        end
        if size(tmp,1) > 1
             call = [call; mean(tmp)];
        end         
    end
end

figure
plot(freq,smoothvect(call(4,:),g),'k')
set(gca,'xlim',[2 139],'ylim',[0.4 0.8],'FontSize',18)
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Coherence','FontSize',18)
title('MEC - CA1','FontSize',18)
box off

%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_mecca1_coherence_example.pdf', m, d, y);
print('-dpdf', savestring)

%% Plot coherence for low and high gamma

% Look at coherence during low and high gamma for CA3-CA1
low = [lookup(25,freq) lookup(55,freq)];
high = [lookup(65,freq) lookup(140,freq)];

call = [];
for an = 1:length(c)
    for d = 1:length(c(an).output)
        tmp = [];
        for e = 1:length(c(an).output{d})
            if ~isempty(c(an).output{d})
                  tmp = [tmp; ...
                    c(an).output{d}(e).coherence([1:noise(1)-1 noise(2)+1:range])];
            end
        end
        if size(tmp,1) > 1
             call = [call; mean(tmp)];
        end         
    end
end

y = [sum(call(:,low(1):low(2)),2)./(low(2)-low(1)) sum(call(:,high(1):high(2)),2)./(high(2)-high(1))]';
figure
plot([1 2],y,'k')
set(gca,'xlim',[0.8 2.2],'ylim',[0.35 0.7],'xTick',[1 2], ...
    'xTickLabel',['Slow Gamma'; 'Fast Gamma'],'FontSize',18)
ylabel('Coherence','FontSize',18)
title({'CA3 - CA1'; '80 out of 85 recording pairs show higher coherence during slow gamma'},'FontSize',18)
box off

%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_ca3ca1_coherence_groupdata.pdf', m, d, y);
print('-dpdf', savestring)

% Look at coherence during low and high gamma for MEC-CA1
call = [];
for an = 1:length(m)
    for d = 1:length(m(an).output)
        tmp = [];
        for e = 1:length(m(an).output{d})
            if ~isempty(m(an).output{d})
                if any(m(an).output{d}(e).index(4) == [6 7])
                 tmp = [tmp; ...
                    m(an).output{d}(e).coherence([1:noise(1)-1 noise(2)+1:range])];
                end
            end
        end
        if size(tmp,1) > 1
             call = [call; mean(tmp)];
        end         
    end
end

y = [sum(call(:,low(1):low(2)),2)./(low(2)-low(1)) sum(call(:,high(1):high(2)),2)./(high(2)-high(1))]';
figure
plot([1 2],y,'k')
set(gca,'xlim',[0.8 2.2],'ylim',[0.35 0.7],'xTick',[1 2], ...
    'xTickLabel',['Slow Gamma'; 'Fast Gamma'],'FontSize',18)
ylabel('Coherence','FontSize',18)
title({'MEC - CA1'; '19 out of 22 recording pairs show higher coherence during slow gamma'},'FontSize',18)
box off

%--------------------------------------------------------------------------
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_mecca1_coherence_groupdata.pdf', m, d, y);
print('-dpdf', savestring)

