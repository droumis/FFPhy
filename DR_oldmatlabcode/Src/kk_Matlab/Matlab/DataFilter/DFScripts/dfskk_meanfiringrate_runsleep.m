
plot_only = 0;

if plot_only~=1
%Animal selection
%-----------------------------------------------------
animals = {'Chapati','Egypt'};

%Filter creation
%--------------------------------------------------------


% epoch filter
epochfilter{1} = ['isequal($type, ''sleep'')'];
epochfilter{2} = ['isequal($type, ''run'')'];


% cell filter
ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''principal''))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''principal''))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''principal''))';

% time filter
%timefilter = { {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 0))', 6} };

% iterator
iterator = 'singlecellanal';

% filter creation
ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter, 'iterator', iterator);
ca2f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca2cellfilter, 'iterator', iterator);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter, 'iterator', iterator);

% set analysis function
ca1f = setfilterfunction(ca1f, 'calctotalmeanrate', {'spikes'});
ca2f = setfilterfunction(ca2f, 'calctotalmeanrate', {'spikes'});
ca3f = setfilterfunction(ca3f, 'calctotalmeanrate', {'spikes'});

% run analysis
ca1f = runfilter(ca1f);
ca2f = runfilter(ca2f);
ca3f = runfilter(ca3f);

end





% Plot

%%% super:   {CA region}{animal}/output/{state}
super{1}=ca1f;
super{2}=ca2f;
super{3}=ca3f;

%%% Matched (throw out units that only exist in either RUN or SLEEP) %%%%%%%

for r=1:3
     % First reconstitute the day-epoch-tetrode-cell index in .index field
 for a=1:length(animals)
    super{r}(a).index=super{r}(a).data;                 % mirror .data field
    for i=1:2                                   % run and sleep
        for j=1:length(super{r}(a).data{i})        % iterate over day-epoch
            dayepoch=super{r}(a).epochs{i}(j,:);
            super{r}(a).index{i}{j}=horzcat(repmat(dayepoch,size(super{r}(a).index{i}{j},1),1), ...
                super{r}(a).index{i}{j});
        end
        dummy=[];
        for j=1:length(super{r}(a).data{i})
            if isempty(super{r}(a).index{i}{j})
                disp('warning means day has no cells..')
            end
            dummy=[dummy ; super{r}(a).index{i}{j}];
        end
        super{r}(a).index{i}=dummy;
    end
 end


% Create .matchedoutput  and .matchedindex :
    % This pairs RUN vs SLEEP epochs & ignores the rest.
for a=1:length(animals)
            super{r}(a).matchedoutput{1} = [];       % run
            super{r}(a).matchedoutput{2} = [];       % sleep
            super{r}(a).matchedindex{1} =  [];       
            super{r}(a).matchedindex{2} =  [];         
            rundummyindex=super{r}(a).index{2};
    for i=1:size(super{r}(a).index{1},1)
        probeindex=super{r}(a).index{1}(i,:);
        rowind = rowfind(probeindex([1 3 4]),rundummyindex(:,[1 3 4]));   % look for matching cell (d-t-c)
        if rowind ~= 0
                super{r}(a).matchedoutput{1} = vertcat(super{r}(a).matchedoutput{1},super{r}(a).output{1}(i));      
                super{r}(a).matchedoutput{2} = vertcat(super{r}(a).matchedoutput{2},super{r}(a).output{2}(rowind)); 
                super{r}(a).matchedindex{1} = vertcat(super{r}(a).matchedindex{1},probeindex);
                super{r}(a).matchedindex{2} = vertcat(super{r}(a).matchedindex{2},super{r}(a).index{2}(rowind,:));
                rundummyindex(rowind,:)=[0 0 0 0];      % discount epoch
        end
    end
end
end








%%% Unmatched, all data %%%

edges=0:0.5:60;
 
% CA1
N=[];
figure
hold on
    %sleep (black)
    dummy=[];
    for i=1:length(animals)
        if ~isempty(super{1}(i).output)
            dummy=[dummy ; super{1}(i).output{1}];
        end
    end
N(1,:) = histc(dummy,edges);
    % run (grey)
    dummy=[];
    for i=1:length(animals)
        if ~isempty(super{1}(i).output)
            dummy=[dummy ; super{1}(i).output{2}];
        end
    end
N(2,:) = histc(dummy,edges);
    b = bar(edges,N','FaceColor',[0.8 0.8 0.8]);
    legend 'RUN';
    axis tight;
    ch=get(b,'children');
    set(ch{1},'FaceColor','k')
    leg = legend([ch{1} ch{2}],'SLEEP','RUN','FontWeight','bold');
    leg = findobj(leg,'type','text');
    set(leg,'FontSize',11,'FontWeight','bold')
    xlabel('mean firing rate (Hz)')
    ylabel('# of units')
    title('CA1','FontSize',30,'FontWeight','bold','FontWeight','bold')
    

% CA2
N=[];
figure
hold on
    %sleep (black)
    dummy=[];
    for i=1:length(animals)
        if ~isempty(super{2}(i).output)
            dummy=[dummy ; super{2}(i).output{1}];
        end
    end
N(1,:) = histc(dummy,edges);
    % run (grey)
    dummy=[];
    for i=1:length(animals)
        if ~isempty(super{2}(i).output)
            dummy=[dummy ; super{2}(i).output{2}];
        end
    end
N(2,:) = histc(dummy,edges);
    b = bar(edges,N','FaceColor',[0.8 0.8 0.8]);
    legend 'RUN';
    axis tight;
    ch=get(b,'children');
    set(ch{1},'FaceColor','k')
    leg = legend([ch{1} ch{2}],'SLEEP','RUN','FontWeight','bold');
    leg = findobj(leg,'type','text');
    set(leg,'FontSize',11,'FontWeight','bold')
    xlabel('mean firing rate (Hz)')
    ylabel('# of units')
    title('CA2','FontSize',30,'FontWeight','bold','FontWeight','bold')
    
    
    
% CA3
N=[];
figure
hold on
    %sleep (black)
    dummy=[];
    for i=1:length(animals)
        dummy=[dummy ; super{3}(i).output{1}];
    end
N(1,:) = histc(dummy,edges);
    % run (grey)
    dummy=[];
    for i=1:length(animals)
        dummy=[dummy ; super{3}(i).output{2}];
    end
N(2,:) = histc(dummy,edges);
    b = bar(edges,N','FaceColor',[0.8 0.8 0.8]);
    legend 'RUN';
    axis tight;
    ch=get(b,'children');
    set(ch{1},'FaceColor','k')
    leg = legend([ch{1} ch{2}],'SLEEP','RUN','FontWeight','bold');
    leg = findobj(leg,'type','text');
    set(leg,'FontSize',11,'FontWeight','bold')
    xlabel('mean firing rate (Hz)')
    ylabel('# of units')
    title('CA3','FontSize',30,'FontWeight','bold','FontWeight','bold')