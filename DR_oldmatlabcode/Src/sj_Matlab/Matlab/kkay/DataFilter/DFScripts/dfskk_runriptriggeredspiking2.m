
% this function can be used for consensus ripples



runscript = 1;

if runscript == 1

% Animal Selection
animals = {'Egypt'};

% Epoch Filter
epochfilter = '(isequal($type, ''run'') || isequal($type, ''sleep''))';

% Time Filter
timefilter = { {'get2dstate', '$velocity < 4'}, ...
               {'kk_getriptimes', '($nripples >= 3)', [], 'cellfilter', '(isequal($area, ''CA1''))'}};

% 
ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''interneuronal''))';    %% ($meanrate < 3))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''interneuronal''))';   
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''interneuronal''))';

%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 3))';
%ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && ($meanrate < 3))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 3))';

% Iterator
iterator = 'kk_multicellanal';

% Filter Creation
ca1f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca1cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca2f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca2cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca3f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca3cellfilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
ca1f = setfilterfunction(ca1f, 'dfakk_getriptriggeredspiking2', {'ripple','ripples','spikes','tetinfo'},'tetfilter','(isequal($area, ''CA1''))', ...
                                                                                        'minthresh',3,'window',[0.5 0.5]);
ca2f = setfilterfunction(ca2f, 'dfakk_getriptriggeredspiking2', {'ripple','ripples','spikes','tetinfo'},'tetfilter','(isequal($area, ''CA1''))', ...
                                                                                        'minthresh',3,'window',[0.5 0.5]);
ca3f = setfilterfunction(ca3f, 'dfakk_getriptriggeredspiking2', {'ripple','ripples','spikes','tetinfo'},'tetfilter','(isequal($area, ''CA1''))', ... 
                                                                                        'minthresh',3,'window',[0.5 0.5]);


% Run Analysis
ca1f = runfilter(ca1f);
ca2f = runfilter(ca2f);
ca3f = runfilter(ca3f);

end





% Plot.  

% CA1-CA1

for reg=[1:3]       % iterate over regions
    if reg==1
        region='CA1';
        clr='k';
        f=ca1f;
    elseif reg==2
        region='CA2';
        clr='g';
        f=ca2f;
    else
        region='CA3';
        clr='r';
        f=ca3f;
    end

% Consolidate cells across epochs, within a day

%collect all cell indices, ignoring specific epochs, into daytetcell
indices = [];
for i=1:length(f.output{1})  % iterate over epochs
    indices = [indices ; f.output{1}(i).cellindices];
end
daytetcell = unique(indices(:,[1 3 4]),'rows');

f.celloutput = struct;

for ind=1:size(daytetcell,1)
    f.celloutput(ind).daytetcell = daytetcell(ind,:);
    f.celloutput(ind).c1vsc2 = [];
    f.celloutput(ind).nospikes = 0;
    f.celloutput(ind).noripples = 0;
    f.celloutput(ind).noepochs = 0;
    %consolidate data across epochs
    for ep=1:length(f.output{1})
        if ~isempty(f.output{1}(ep).cellindices)   % epoch may not have any data (fully excluded)
            jj=rowfind(daytetcell(ind,:),f.output{1}(ep).cellindices(:,[1 3 4]));
            if jj
                f.celloutput(ind).c1vsc2 = [f.celloutput(ind).c1vsc2 ; f.output{1}(ep).c1vsc2(jj,:)];
                f.celloutput(ind).nospikes = f.celloutput(ind).nospikes + f.output{1}(ep).nospikes(jj);
                f.celloutput(ind).noripples = f.celloutput(ind).noripples + f.output{1}(ep).noripples(jj);
                f.celloutput(ind).noepochs = f.celloutput(ind).noepochs + 1;
                f.celloutput(ind).time = f.output{1}(ep).time;   % (all time vectors are the same..)
            end
        end
    end
    f.celloutput(ind).c1vsc2total = sum(f.celloutput(ind).c1vsc2,1);
    f.celloutput(ind).c1vsc2totalnorm = f.celloutput(ind).c1vsc2total/sqrt(f.celloutput(ind).nospikes*f.celloutput(ind).noripples);
end
    
    
    % plot unnormalized unit cross-correlograms
    if 0
    for c=1:length(f.celloutput)
        if counter > 60
            figure
            counter = 1;
            titlestring=sprintf('%s %s ripple triggered spiking',animals{1},region);
            title(titlestring)
        end
        subplot(6,10,counter)
            h = bar(f.celloutput(c).time,sum(f.celloutput(c).c1vsc2(1,:),1));
            set(h(1),'facecolor',[1 1 1])  
        titlestring=sprintf('%d %d %d',f.celloutput(c).daytetcell);
        title(titlestring)
        
        counter = counter+1;  
    end
    end
    
    
    
    super{reg} = f;
   
    
end



% plot regional aggregates (check if normalized or not)

        figure
        hold on
clear c1vsc2_reg

if 1
    
for reg=[1 2 3]
    time = super{reg}.celloutput(end).time;
    c1vsc2_reg{reg} = [];
    for c=1:length(super{reg}.celloutput)
        %c1vsc2_reg{reg} = [c1vsc2_reg{reg} ; super{reg}.celloutput(c).c1vsc2total];
        c1vsc2_reg{reg} = [c1vsc2_reg{reg} ; super{reg}.celloutput(c).c1vsc2total];
    end
    h = bar(time,sum(c1vsc2_reg{reg},1));
    if reg==1
        set(h(1),'facecolor',[0 0 1])
        alpha(0.5)
    elseif reg==3
        set(h(1),'facecolor',[1 0 0])
        alpha(0.5)
    else
        set(h(1),'facecolor',[0 1 0])
        alpha(0.5)        
    end
end

end

% plot regional aggregates (normalized)












