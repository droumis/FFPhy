%% Define Filter Options

% Animal selection
animals = {'Bond','Conley','Corriander','Dudley','Eight','Five','Frank','Miles','Ten'};
%animals = {'Dudley','Miles','Bond','Frank','Ten'};
%animals = {'Corriander','Frank'};

% Epoch selection
epochfilter = [];
for i = 1:24
    epochfilter{i} = ['($experimentday == ',num2str(i),')'];
end

% Tetrode selection
ca1tetfilter = '(isequal($area, ''CA1'') & ($numcells > 2))';

% Time selection
timefilter = {{'getlinstate', 'abs($velocity) > 3', 1}};

%Select iterator
iterator = 'eeganal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter, 'eegtetrodes',ca1tetfilter,'iterator', iterator);

f = setfilterfunction(f,'modulationindex',{'eeg'});
f = runfilter(f);

%% Plot MI
center = (-1*pi+pi/20):pi/10:(pi-pi/20);
for an = 1:length(f)
    for d=1:length(f(an).output)
        for e=1:length(f(an).output{d})
            %figure
            imagesc([center 2*pi+center],20:5:200,[f(an).output{d}(e).modulationindex ...
            f(an).output{d}(e).modulationindex]);
            axis xy
            set(gca,'clim',[0.045 0.055])
            f(an).output{d}(e).index
            pause
         end
        
    end
       
end

%% PLOT MI_2
for an = 1:length(f)
    temp = zeros(61,length(f(an).output));
    for d=1:length(f(an).output)
        for e=1:length(f(an).output{d})
            temp(:,d) = [temp(:,d)+f(an).output{d}(e).mi];
        end
        temp(:,d) = temp(:,d)/length(f(an).output{d});
    end  
    figure(an)
    imagesc([],20:5:200,temp);
    axis xy
end

%% NEW ATTEMPTS
for an=1:length(f)
    temp = [];
    count = 0;
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            count = count + 1;
            if isempty(temp)
                temp = f(an).output{d}(e).modulationindex;
            else
                temp = temp+f(an).output{d}(e).modulationindex;
            end
        end
    end
    temp = temp./count;
    figure
    imagesc([center 2*pi+center],20:5:200,[temp temp])
    axis xy
    %set(gca,'clim',[0.0475 0.0525])
end

g = 20:2:140;
for an=1:length(f)
    tmp = [];
    count = 0;
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            count = count + 1;
            if isempty(tmp)
                tmp = f(an).output{d}(e).mi;
            else
                tmp = tmp+f(an).output{d}(e).mi;
            end
            count = count+1;
        end
    end
    [maxtab mintab] = peakdet(tmp,0.00000001);
    figure(an)
    plot(g,tmp,g(mintab(:,1)),mintab(:,2),'g*')
end

for an=1:length(f)
    for d = 1:length(f(an).output)
        if ~isempty(f(an).output{d})
            temp = [];
            count = 1;
            for e = 1:length(f(an).output{d})
                if isempty(temp)
                    temp = f(an).output{d}(e).modulationindex;
                else
                    temp = temp+f(an).output{d}(e).modulationindex;
                    count = count +1;
                end
            end
            temp = temp./count;
            imagesc([center 2*pi+center],20:5:200,[temp temp])
            axis xy
            f(an).animal(1), f(an).epochs{d}
            pause
        end
    end
end
