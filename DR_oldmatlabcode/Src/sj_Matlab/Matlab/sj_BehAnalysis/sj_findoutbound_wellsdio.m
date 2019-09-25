
function [outbound_curve, all_outbound_logic, all_outbound_wellstend, ntrajs_perday] = sj_findoutbound_wellsdio (days,epochs,prefix, doprobfit)

%% sj_findoutbound_wellsdio ([1:3,6:7],[2 4], 'RCa', 1);
%% Calls sj_day_findinbound by looping over all days, and returns vector of
%% outbound_logic (Correct incorrect), wells exited and entered, and number
%% of outbound trajectories per day
%% Shantanu Jadhav, 12/30/09

if nargin<=3,
    doprobfit = 0;
end

%%Initialize
all_outbound_logic = [];
all_outbound_wellstend = [];
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
        
        %% Find outbound
        outbound_stidx = find(wells(:,1)==1);
        outbound_wellstend = wells(outbound_stidx,:);
        
        %% Find out which outbound is correct, and return the logic
        outbound_logic = zeros(length(outbound_wellstend),1);
        mintr = min(outbound_stidx);
        if mintr==1,
            outbound_stidx(find(outbound_stidx==1))=[];
        end
        
        %% Find outbound
        outbound_stidx = find(wells(:,1)==1);
        outbound_wellstend = wells(outbound_stidx,:);
        
        %% Find out which outbound is correct, and return the logic
        outbound_logic = zeros(length(outbound_wellstend),1);
        mintr = min(outbound_stidx);
        if mintr==1,
            outbound_stidx(find(outbound_stidx==1))=[];
        end
        
        corridx = find( (wells(outbound_stidx,2)~=wells(outbound_stidx-1,1))...
            & (wells(outbound_stidx,2)~=1) );
        correct = outbound_stidx(corridx);
        
        
        %% Taking into account trials in which he goes from well 1 to well 1
        
        % What is this? dont bother
%         retrace = correct(find (wells(correct-1,1)==1));
%         retrace_idx = find (wells(correct-1,1)==1);
%         
%         % Only bother if its atleast 3 trajectories deep
%         remove_retr=[];
%         for i=1:length(retrace)
%             if (retrace(i)>=3) & (wells(retrace(i),2) == wells(retrace(i)-2,1))
%                 remove_retr = [remove_retr; retrace_idx(i)];
%                 %remove_retr = retrace(find(wells(retrace(i),2)==wells(retrace(i)-2,1)));
%                 %corr(remove_retr)=[];
%             end
%         end
%         correct(remove_retr)=[];
        
        %% Either get rid of first trial, or keep it in as correct
        % Alt 1 - Keep it as correct
        % if mintr==1.
        %     correct = [1;correct]; outbound_stidx = [1;outbound_stidx];
        % end
        % Alt 2 - Remove it - Do nothing and update length of outbound_logic
        outbound_logic = zeros(length(outbound_wellstend)-1,1);
        
        [c,ia,ib] = intersect(outbound_stidx, correct);
        
        outbound_logic(ia) = 1;
        sum(outbound_logic);
        
        
        %% Add to cumulative count
        all_outbound_logic = [all_outbound_logic; outbound_logic];
        all_outbound_wellstend = [all_outbound_wellstend; outbound_wellstend];
        n_temp = n_temp + length(outbound_logic);
    end
    cnt=cnt+1;
    ntrajs_perday(cnt) = n_temp;
    
end

%savefile = [prefix '_inbound'];
%save(savefile,'all_inbound_logic', 'all_inbound_wellstend', 'ntrajs_perday');


%%% Do Algorithm

if doprobfit == 1
    
    [pc, lt] = getestprobcorrect(all_outbound_logic, 0.5, 0);
    outbound_curve = pc(2:end,:);
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
    title([prefix ' - Outbound Trials'],'FontSize',24,'Fontweight','bold');
    
    %% Axes Names
    xlabel('Trial Number')
    ylabel('Probability of a Correct Response')
    
end
%savefile = [prefix '_outbound'];
%save(savefile,'all_outbound_logic',
%'all_outbound_wellstend','ntrajs_perday','outbound_curve');