
function [inbound_curve, all_inbound_logic, all_inbound_wellstend, ntrajs_perday] = sj_findinbound_wellsdio (days,epochs,prefix, doprobfit)

%% sj_findinbound_wellsdio ([1:4,6:7],[2 4], 'RCa', 1);
%% Calls sj_day_findinbound by looping over all days, and returns vector of
%% inbound_logic (Correct incorrect), wells exited and entered, and number
%% of inbound trajectories per day
%% Shantanu Jadhav, 12/30/09

if nargin<=3,
    doprobfit = 0;
end

%%Initialize
all_inbound_logic = [];
all_inbound_wellstend = [];
ntrajs_perday = [];

%% Loop over days
cnt=0;
for day=days
    
    dsz = '';
    if (day < 10)
       dsz = '0';
    end
    
    filename = [prefix 'wellsdio' dsz num2str(day)];
    load(filename);
    
    n_temp=0;
    for epoch = epochs
        
        wells=wellsdio{day}{epoch};
        %% Find inbound
        inbound_stidx = find(wells(:,1)~=1);
        inbound_wellstend = wells(inbound_stidx,:);
        %% Find out which inbound is correct, and return the logic
        inbound_logic = zeros(length(inbound_wellstend),1);
        corr = find(inbound_wellstend(:,2)==1);
        inbound_logic(corr) = 1;
        
        %[inbound_logic, inbound_wellstend] = sj_day_findinbound (linpos,i,epoch);
        
        %% Add to cumulative count
        all_inbound_logic = [all_inbound_logic; inbound_logic];
        all_inbound_wellstend = [all_inbound_wellstend; inbound_wellstend];
        n_temp = n_temp + length(inbound_logic);
    end
    cnt=cnt+1;
    ntrajs_perday(cnt) = n_temp;
    
end


%%% Do Algorithm

if doprobfit == 1
    
    [pc, lt] = getestprobcorrect(all_inbound_logic, 0.5, 0);
    inbound_curve = pc(2:end,:);
    figure; hold on;
    
    %redimscreen;
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
    set(0,'defaultaxeslinewidth',2);
    
    %% Get time (day) axis
    t = 1:size(pc,1)-1;
    %%% Plot data %%
    plot(t', pc(2:end,1),'r-', 'Linewidth', 2);
    plot(t', pc(2:end,2),'m--', 'Linewidth', 2);
    plot(t', pc(2:end,3),'m--', 'Linewidth', 2);
    
    %% Plot trials in a day
    shift = [0, ntrajs_perday(1:end-1)];
    cumshift = cumsum(shift);
    
    for i=1:length(days)
        %line([ntrajs_perday(i) ntrajs_perday(i)], [0 1],'g--');
        plot((ntrajs_perday(i)+cumshift(i))*ones(size([0:0.1:1])) , [0:0.1:1], 'k--');
    end
    
    %% Axis
    axis([1 t(end)  0 1]);
    
    %% Backgnd Prob
    %line([1 t(end)], [background_prob  background_prob ], 'k');
    plot(t, 0.5*ones(size(t)),'k', 'Linewidth', 2);
    
    %% Title
    title([prefix ' - Inbound Trials'],'FontSize',24,'Fontweight','bold');
    
    %% Axes Names
    xlabel('Trial Number')
    ylabel('Probability of a Correct Response')
    
end

%savefile = [prefix '_inbound'];
%save(savefile,'all_inbound_logic',
%'all_inbound_wellstend','ntrajs_perday','inbound_curve');