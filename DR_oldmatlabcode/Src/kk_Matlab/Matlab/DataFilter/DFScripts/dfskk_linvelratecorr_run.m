
runscript = 1;


if runscript==1

for s=1:3
    
%Animal selection
%-----------------------------------------------------
animals = {'Chapati'};

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
   str = '(($traj ~= -1) & (abs($velocity) <= 0.5))';
elseif s==2
   str = '(($traj ~= -1) & (abs($velocity) >= 8))';
else
   str = '(($traj ~= -1) & (abs($velocity) > 0.5) & (abs($velocity) < 8))';
end
    
timefilter = { {'getlinstate', str, 6} };

% iterator
iterator = 'singlecellanal';

% filter creation
ca1f{s} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime',timefilter,'iterator',iterator);
ca2f{s} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca2cellfilter,'excludetime',timefilter,'iterator',iterator);
ca3f{s} = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime',timefilter,'iterator',iterator);

% set analysis function
ca1f{s} = setfilterfunction(ca1f{s}, 'calclinvelratecorr', {'spikes','linpos'});
ca2f{s} = setfilterfunction(ca2f{s}, 'calclinvelratecorr', {'spikes','linpos'});
ca3f{s} = setfilterfunction(ca3f{s}, 'calclinvelratecorr', {'spikes','linpos'});

% run analysis

ca1f{s} = runfilter(ca1f{s});
ca2f{s} = runfilter(ca2f{s});
ca3f{s} = runfilter(ca3f{s});

end

end




% Plot
edges=0:0.1:1;

% CA1
N=[];
figure
hold on
N(1,:) = histc(ca1f{1}.output{1},edges);   % state1 (black)
N(2,:) = histc(ca1f{2}.output{1},edges);   % state2 (white)
    b = bar(edges,N','FaceColor','k');
    legend 'State1' 'State2';
    axis tight;
    ch=get(b,'children');
    set(ch{1},'FaceColor','k')
    set(ch{2},'FaceColor','w','EdgeColor','k')
    leg = legend([ch{1} ch{2}],'state1','state2','FontWeight','bold');
    leg = findobj(leg,'type','text');
    set(leg,'FontSize',11,'FontWeight','bold')
    xlabel('mean firing rate (Hz)')
    ylabel('# of units')
    title('CA1','FontSize',30,'FontWeight','bold')
    subtitle = {'peak firing','run epochs','state (lin)'};
    gtext(subtitle,'HorizontalAlignment','center','FontSize',15)
    
% CA2
N=[];
figure
hold on
N(1,:) = histc(ca2f{1}.output{1},edges);   % state1 (black)
N(2,:) = histc(ca2f{2}.output{1},edges);   % state2 (white)
    b = bar(edges,N','FaceColor','k');
    legend 'State1' 'State2';
    axis tight;
    ch=get(b,'children');
    set(ch{1},'FaceColor','k')
    set(ch{2},'FaceColor','w','EdgeColor','k')
    leg = legend([ch{1} ch{2}],'state1','state2','FontWeight','bold');
    leg = findobj(leg,'type','text');
    set(leg,'FontSize',11,'FontWeight','bold')
    xlabel('mean firing rate (Hz)')
    ylabel('# of units')
    title('CA2','FontSize',30,'FontWeight','bold')
    subtitle = {'peak firing','run epochs','state (lin)'};
    gtext(subtitle,'HorizontalAlignment','center','FontSize',15)
    
% CA3
N=[];
figure
hold on
N(1,:) = histc(ca3f{1}.output{1},edges);   % state1 (black)
N(2,:) = histc(ca3f{2}.output{1},edges);   % state2 (white)
    b = bar(edges,N','FaceColor','k');
    legend 'State1' 'State2';
    axis tight;
    ch=get(b,'children');
    set(ch{1},'FaceColor','k')
    set(ch{2},'FaceColor','w','EdgeColor','k')
    leg = legend([ch{1} ch{2}],'state1','state2','FontWeight','bold');
    leg = findobj(leg,'type','text');
    set(leg,'FontSize',11,'FontWeight','bold')
    xlabel('mean firing rate (Hz)')
    ylabel('# of units')
    title('CA3','FontSize',30,'FontWeight','bold')
    subtitle = {'peak firing','run epochs','state (lin)'};
    gtext(subtitle,'HorizontalAlignment','center','FontSize',15)
    %h = text(5,5,subtitle,'HorizontalAlignment','center');
    %    set(h,'FontSize',15)
    