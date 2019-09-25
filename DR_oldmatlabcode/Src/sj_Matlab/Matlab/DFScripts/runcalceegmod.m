% Animal selection
animals = {'Bond','Conley','Corriander','Dudley','Eight','Five','Frank','Miles','Six','Ten'};

% Epoch selection
epochfilter = [];
for i = 1:24
    epochfilter{i} = ['($experimentday == ',num2str(i),')'];
end

% Tetrode selection
ca1tetfilter = '(isequal($area, ''CA1'')) & $numcells>=1';

% Cell Selection
ca1cellfilter = '(isequal($area, ''CA1'') & ($numspikes > 200))';

% Time selection
timefilter = {{'getgammatimes', '($ngamma == 1)', [], 'high','tetfilter', ca1tetfilter}, ...
    {'getlinstate', 'abs($velocity) >3', 6}};

% Iterator selection
iterator = 'singlecelleeganal';


% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime', ...
    timefilter,'cells', ca1cellfilter, 'eegtetrodes',ca1tetfilter,'iterator', iterator);

f = setfilterfunction(f, 'calceegmod', {'highgamma','spikes'},'appendindex',1,'nbins',30);
%f = setfilterfunction(f, 'calcthetamod', {'theta','spikes'},'appendindex',1,'nbins',30);

f = runfilter(f);

%% Plot

%Try plotting each cell
bins = -pi:2*pi/20:pi;
bin = [bins(1:end-1) bins(1:end-1)+2*pi];
R = cell(24,1);
A = cell(24,1);
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if length(f(an).output{d}(e).phase) > 100
                temp = f(an).output{d}(e).phase;    
                day = f(an).output{d}(e).index(1);
                % To compute the mean direction a 
                x = sum(cos(temp));
                y = sum(sin(temp));

                if x > 0 && y >0
                    a = atan(y/x);
                elseif x < 0
                    a = atan(y/x)+pi;
                elseif x> 0 && y < 0
                    a = atan(y/x)+2*pi;
                end
                A{day} = [A{day} a];
                % To compute mean length
                r = sqrt(x.^2 + y.^2)./length(temp);
                R{day} = [R{day} r];
%                 count = histc(f(an).output{d}(e).phase,bins)./length(f(an).output{d}(e).phase);
%                 bar(bin,[count(1:end-1); count(1:end-1)])
%                 set(gca,'ylim',[0 0.1])
%                title(num2str([an f(an).output{d}(e).index]))
%                pause
            end
        end
    end
end

r = [];
a = [];
for d = 1:length(R)
    r = [r; d*ones(length(R{d}),1) R{d}'];
    a = [a; d*ones(length(A{d}),1) A{d}'];
end
plot(r(:,1),r(:,2),'.',a(:,1),a(:,2),'r.')
[b,bint,r,rint,stats] = regress(r(:,2),[ones(length(r),1) r(:,1)],0.05);


%Plot each epoch separately
bins = -pi:(2*pi/20):pi;
for an = 1:length(f)
    tmp = structure_group_combine(f(an),'output','phase');
    for d = 1:length(tmp)
        if length(tmp{d})>500
            count = histc(tmp{d}, bins)./length(tmp{d});
            figure(an)
            hold on
            plot([bins(1:end-1) bins(1:end-1)+2*pi],[count(1:end-1); count(1:end-1)])
            set(gca,'ylim',[0 0.1])
            title(num2str([an d]))
        end
    end
end


%Combine all epochs
bins = -pi:(2*pi/11):pi;
for an = 1:length(f)
    temp = [];
    tmp = structure_group_combine(f(an),'output','phase');
    for d = 1:length(tmp)
        temp = [temp; tmp{d}];
    end
    count = histc(temp, bins)./length(temp);
    figure(an)
    bar([bins(1:end-1) bins(1:end-1)+2*pi],[count(1:end-1); count(1:end-1)])
    set(gca,'ylim',[0 0.15])
end


%To look at mean direction over time
figure
hold on
for an=1:length(f)
    tmp = structure_group_combine(f(an),'output','phase');
    R = nan(length(tmp),1);
    V = nan(length(tmp),1);
    for d = 1:length(tmp)
        if length(tmp{d})>500
            temp = tmp{d};    

            % To compute the mean direction a 
            x = sum(cos(temp));
            y = sum(sin(temp));

            if x > 0 && y >0
                a = atan(y/x);
            elseif x < 0
                a = atan(y/x)+pi;
            elseif x> 0 && y < 0
                a = atan(y/x)+2*pi;
            end

            % To compute mean length
            r = sqrt(x.^2 + y.^2)./length(temp);
            R(d) = r;
            % To compute circular standard deviation
            v =sqrt(-2*log(r));
            V(d) = v;
        end
    end
plot(1:d,R,'b*',1:d,V,'r*')
end
