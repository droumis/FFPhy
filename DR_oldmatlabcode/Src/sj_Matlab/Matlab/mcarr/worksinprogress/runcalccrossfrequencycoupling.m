% Animal selection
animals = {'Bond','Frank','Ten'};

% Epoch selection
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
%epochfilter{2} = 'isequal($type,''sleep'')';

% Tetrode pair selection
tetfilter = {'isequal($area, ''CA1'') & $maxcell==1','isequal($area, ''CA3'')&$maxcell==1'};

%Iterator Selection
iterator = 'epocheegnonreferenceanal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter, 'eegtetrodepairs',tetfilter,'iterator', iterator);
f = setfilterfunction(f, 'calccrossfrequencycoupling', {'ripples','cellinfo','ripple','lowgamma'});
f = runfilter(f);

%% PLOT

% Look using HR method:
    % Lakatos, P. et al. (2005) An oscillatory hierarchy controlling 
    % neuronal excitability and stimulus processing in the auditory cortex.
    % J. Neurophysiology 94, 1904–1911
    
% First plot an example demonstrating relationship between gamma phase and
% ripple amplitude. Three cycles shown for illustration. Ripple amplitude
% in CA1 is modulated by gamma oscillations in CA3. Same epoch as other
% examples (2a,2b,3a,and 3d)

bin = -pi:pi/6:pi-pi/6;
bin = bin*10000;
hr = accumarray(lookup(f(3).output{1}(9).phase,bin),f(3).output{1}(9).amp,[length(bin) 1],@(x) mean(x));
hr_err = accumarray(lookup(f(3).output{1}(9).phase,bin),f(3).output{1}(9).amp,[length(bin) 1],@(x) stderr(x(:)));
gam = accumarray(lookup(f(3).output{1}(9).phase,bin),f(3).output{1}(9).gamma,[length(bin) 1],@(x) mean(x));
U_gam = hr' + hr_err'; U_gam = [U_gam U_gam];
L_gam = hr'-hr_err';  L_gam = [L_gam L_gam];

figure
bin = -pi:pi/6:pi-pi/6;
bin = [bin 2*pi+bin];
plot(bin,[hr; hr],'k')
hold on
fill([bin fliplr(bin)], [U_gam fliplr(L_gam)],'k','FaceAlpha',1,'EdgeColor','none')
set(gca,'xtick',[-pi 0 pi 2*pi 3*pi],'xticklabel',[180 0 180 0 180 0 180],'ylim',[1 1.5])
xlabel('CA3 Gamma phase')
ylabel('CA1 Ripple amplitude')
legend([{'CA1 Ripple amplitude'},{'CA3 Average Gamma Trace'}])

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_cross_frequency_coupling_example.pdf', m, d, y);
print('-dpdf', savestring)


% Determine phase of maximal ripple amplitude for each session
bin = -pi:pi/6:pi-pi/6;
bin = bin*10000;
HR = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).phase)
                if length(f(an).output{d}(e).phase)>10000
                    a = accumarray(lookup(f(an).output{d}(e).phase,bin),f(an).output{d}(e).amp,[length(bin) 1],@(x) mean(x));
                    HR = [HR a];
                end
            end
        end
    end
end
bin = bin/10000;
[maxval maxx] = max(HR); [minval minx] = min(HR);
%Use rayleigh_test to figure out preferred angle
rayleigh_test(bin(maxx))
%theta = 0.8 radians ~= 45 degrees
%p<1e-5
rayleigh_test(bin(minx))
%theta = -2.5 radians = -145 degrees
%p<1e-5

figure
plot(bin,hist(bin(maxx),bin)./length(maxx));
hold on

boxplot(bin(maxx),'orientation','horizontal','color','k','notch','on')
xlabel('CA3 Gamma Phase')
ylabel('Proportion of maximal CA1 ripple amplitude')
set(gca,'xlim',[-pi pi],'ylim',[0 1.25],'ytick',0:0.025:0.25,'yticklabel',0:0.025:0.25,'xtick',bin,'xticklabel',round(bin*180/pi))

%Depth of modulation: average depth of modulation = 15%
q = ((maxval-minval)./(maxval+minval));

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_distribution_crossfrequency_coupling.pdf', m, d, y);
print('-dpdf', savestring)

%% RUN WITH CA1 GAMMA and CA1 RIPPLE AMPLITUDE

% Animal selection
animals = {'Bond','Frank','Ten'};

% Epoch selection
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
%epochfilter{2} = 'isequal($type,''sleep'')';

% Tetrode pair selection
tetfilter = 'isequal($area, ''CA1'') & $maxcell==1';

%Iterator Selection
iterator = 'epocheegnonreferenceanal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter, 'eegtetrodes',tetfilter,'iterator', iterator);
f = setfilterfunction(f, 'calccrossfrequencycoupling', {'ripples','cellinfo','ripple','lowgamma'});
f = runfilter(f);

%% PLOT

bin = -pi:pi/6:pi-pi/6;
bin = bin*10000;
hr = accumarray(lookup(f(3).output{1}(9).phase,bin),f(3).output{1}(9).amp,[length(bin) 1],@(x) mean(x));
hr_err = accumarray(lookup(f(3).output{1}(9).phase,bin),f(3).output{1}(9).amp,[length(bin) 1],@(x) stderr(x(:)));
gam = accumarray(lookup(f(3).output{1}(9).phase,bin),f(3).output{1}(9).gamma,[length(bin) 1],@(x) mean(x));
U_gam = hr' + hr_err'; U_gam = [U_gam U_gam U_gam];
L_gam = hr'-hr_err';  L_gam = [L_gam L_gam L_gam];

figure
bin = -pi:pi/6:pi-pi/6;
bin = [-2*pi+bin bin 2*pi+bin];
plot(bin,[hr; hr; hr],'k',bin,[gam/max(gam); gam/max(gam); gam/max(gam)],'r')
hold on
fill([bin fliplr(bin)], [U_gam fliplr(L_gam)],'k','FaceAlpha',1,'EdgeColor','none')
set(gca,'xtick',[-3*pi -2*pi -pi 0 pi 2*pi 3*pi],'xticklabel',[180 0 180 0 180 0 180],'ylim',[-1 2])
xlabel('CA1 Gamma phase')
ylabel('CA1 Ripple amplitude')
legend([{'CA1 Ripple amplitude'},{'CA1 Average Gamma Trace'}])

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_CA1_cross_frequency_coupling_example.pdf', m, d, y);
print('-dpdf', savestring)


% Determine phase of maximal ripple amplitude for each session
bin = -pi:pi/6:pi-pi/6;
bin = bin*10000;
HR = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).phase)
                if length(f(an).output{d}(e).phase)>10000
                    a = accumarray(lookup(f(an).output{d}(e).phase,bin),f(an).output{d}(e).amp,[length(bin) 1],@(x) mean(x));
                    HR = [HR a];
                end
            end
        end
    end
end
bin = bin/10000;
[maxval maxx] = max(HR); [minval minx] = min(HR);
%Use rayleigh_test to figure out preferred angle
rayleigh_test(bin(maxx))
%theta = 0.88 radians ~= 50 degrees
%p<1e-5
rayleigh_test(bin(minx))
%theta = -2.01 radians = -133 degrees
%p<1e-5

figure
plot(bin,hist(bin(maxx),bin)./length(maxx));
xlabel('CA1 Gamma Phase')
ylabel('Proportion of maximal CA1 ripple amplitude')
set(gca,'xlim',[-pi pi],'ylim',[0 0.25],'ytick',0:0.025:0.25)

%Depth of modulation: average depth of modulation = 15%
q = ((maxval-minval)./(maxval+minval));

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_CA1_distribution_crossfrequency_coupling.pdf', m, d, y);
print('-dpdf', savestring)

