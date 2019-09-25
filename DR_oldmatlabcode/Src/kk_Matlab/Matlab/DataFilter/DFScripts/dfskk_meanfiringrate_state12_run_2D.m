
plot_only=1;


if plot_only~=1

for s=1:2
    
%Animal selection
%-----------------------------------------------------
animals = {'Chapati','Egypt'};

%Filter creation
%--------------------------------------------------------

% epoch filter
epochfilter{1} = ['isequal($type, ''run'')'];

% cell filter
ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 100))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($meanrate < 100))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 100))';

% time filter
if s==1
   str = '((abs($velocity) <= 0.1))';
elseif s==2
   str = '((abs($velocity) >= 8))';
else
   str = '(((abs($velocity) > 0.5) & (abs($velocity) < 8))';
end
    
timefilter = { {'get2dstate', str} };

% iterator
iterator = 'singlecellanal';

% filter creation
ca1f{s} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime',timefilter,'iterator',iterator);
ca2f{s} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca2cellfilter,'excludetime',timefilter,'iterator',iterator);
ca3f{s} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime',timefilter,'iterator',iterator);

% set analysis function
ca1f{s} = setfilterfunction(ca1f{s}, 'calctotalmeanrate', {'spikes'});
ca2f{s} = setfilterfunction(ca2f{s}, 'calctotalmeanrate', {'spikes'});
ca3f{s} = setfilterfunction(ca3f{s}, 'calctotalmeanrate', {'spikes'});

% run analysis

ca1f{s} = runfilter(ca1f{s});
ca2f{s} = runfilter(ca2f{s});
ca3f{s} = runfilter(ca3f{s});

end

end


% Plot

edges=0:0.25:10;

% CA1
N=[];
figure
hold on
    %state 1 (black)
    dummy=[];
    for i=1:length(animals)
        dummy=[dummy ; ca1f{1}(i).output{1}];
    end
N(1,:) = histc(dummy,edges);
    % state2 (grey)
    dummy=[];
    for i=1:length(animals)
        dummy=[dummy ; ca1f{2}(i).output{1}];
    end
N(2,:) = histc(dummy,edges);
    b = bar(edges,N','FaceColor',[0.8 0.8 0.8]);
    legend 'RUN';
    axis tight;
    ch=get(b,'children');
    set(ch{1},'FaceColor','k')
    leg = legend([ch{1} ch{2}],'State1','State2','FontWeight','bold');
    leg = findobj(leg,'type','text');
    set(leg,'FontSize',11,'FontWeight','bold')
    xlabel('mean firing rate (Hz)')
    ylabel('# of units')
    title('CA1','FontSize',30,'FontWeight','bold','FontWeight','bold')
    
    
    % CA2
N=[];
figure
hold on
    %state 1 (black)
    dummy=[];
    for i=1:length(animals)
        dummy=[dummy ; ca1f{1}(i).output{1}];
    end
N(1,:) = histc(dummy,edges);
    % state2 (grey)
    dummy=[];
    for i=1:length(animals)
        dummy=[dummy ; ca1f{2}(i).output{1}];
    end
N(2,:) = histc(dummy,edges);
    b = bar(edges,N','FaceColor',[0.8 0.8 0.8]);
    legend 'RUN';
    axis tight;
    ch=get(b,'children');
    set(ch{1},'FaceColor','k')
    leg = legend([ch{1} ch{2}],'State1','State2','FontWeight','bold');
    leg = findobj(leg,'type','text');
    set(leg,'FontSize',11,'FontWeight','bold')
    xlabel('mean firing rate (Hz)')
    ylabel('# of units')
    title('CA2','FontSize',30,'FontWeight','bold','FontWeight','bold')
    
    
    
    
    
    
    
    
    % CA3
N=[];
figure
hold on
    %state 1 (black)
    dummy=[];
    for i=1:length(animals)
        dummy=[dummy ; ca1f{1}(i).output{1}];
    end
N(1,:) = histc(dummy,edges);
    % state2 (grey)
    dummy=[];
    for i=1:length(animals)
        dummy=[dummy ; ca1f{2}(i).output{1}];
    end
N(2,:) = histc(dummy,edges);
    b = bar(edges,N','FaceColor',[0.8 0.8 0.8]);
    legend 'RUN';
    axis tight;
    ch=get(b,'children');
    set(ch{1},'FaceColor','k')
    leg = legend([ch{1} ch{2}],'State1','State2','FontWeight','bold');
    leg = findobj(leg,'type','text');
    set(leg,'FontSize',11,'FontWeight','bold')
    xlabel('mean firing rate (Hz)')
    ylabel('# of units')
    title('CA3','FontSize',30,'FontWeight','bold','FontWeight','bold')
    