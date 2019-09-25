
function [perf, smooth_perf, perfvsday] = sj_smoothperf1_in (dir, days,epochs,prefix,days_ex, win, figopt, savefig)
%function [perf, smooth_perf, perfvsday] = sj_smoothperf1(prefix, win, logic, ntrajs_perday, type, figopt)
% eg [perf, smooth_perf, perfvsday] = sj_smoothperf1_out('/data25/sjadhav/RippleInterruption/REd_direct',[1:8],[2 4],'REd',[], 10, 1,0);

if nargin<6,
    win = 10;
end

if nargin<7,
    figopt=0;
end

if nargin<8,
    savefig=0;
end

% ------------------------------
% Figure and Font Sizes
clr = 'r'; % Main line r/b/k
%clr2 = 'b'; % Confidence bounds m/c/k

forppr = 1; 
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/Behavior/IndividualAnimals/RawBehavior/';
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    titlefont = 18;
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    titlefont = 28;
    xfont = 20;
    yfont = 20;
end

% ---------------------------------------

cd(dir);
currdir = pwd;


epflag=0;  % If epoch is not specified, get run epochs from task file
%%Initialize
all_inbound_logic = [];
all_inbound_wellstend = [];
ntrajs_perday = [];
ntrajs_persess = [];

%% Loop over days
cnt=0; cntsess=0;
for day=days
    
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    filename = [prefix 'linpos' dsz num2str(day)];
    load(filename);
    
     % Get Epochs if not specified
    if isempty(epochs)
        epflag=1;
        epochs=[];
        taskfile = [prefix 'task' dsz num2str(day)];
        load(taskfile);
        toteps = length(linpos{day});
        for i = 1:toteps,
            if ~isempty(linpos{day}{i})
                if ( strcmp('TrackA',task{day}{i}.description) | strcmp('Track A',task{day}{i}.description) )              
                    epochs=[epochs, i];
                end
            end
        end
    end 
    
    
    n_temp=0;
    for e = 1:length(epochs)
        epoch = epochs(e);
        [inbound_logic, inbound_wellstend] = sj_day_findinbound (linpos,day,epoch);
        
        %% Add to cumulative count
        all_inbound_logic = [all_inbound_logic; inbound_logic];
        all_inbound_wellstend = [all_inbound_wellstend; inbound_wellstend];
        
        cntsess=cntsess+1;
        ntrajs_persess(cntsess) = length(inbound_logic);
        
        n_temp = n_temp + length(inbound_logic);
    end
    cnt=cnt+1;
    ntrajs_perday(cnt) = n_temp;
    
     % Reset epochs at end of day
    if epflag==1
        epochs = [];
    end
end

ex_ntrajs_perday = [];
ex_ntrajs_persess = [];
if ~isempty(days_ex)
    ex_inbound_logic = [];
    ex_inbound_wellstend = [];
    
    
    cnt=0; cntsess=0;
    for day=days_ex
        dsz = '';
        if (day < 10)
            dsz = '0';
        end
        
        filename = [prefix 'linpos' dsz num2str(day)];
        load(filename);
        
        n_temp=0;
        for e = 1:length(epochs)
            epoch = epochs(e);
            [inbound_logic, inbound_wellstend] = sj_day_findinbound (linpos,day,epoch);
            
            %% Add to cumulative count
            ex_inbound_logic = [ex_inbound_logic; inbound_logic];
            ex_inbound_wellstend = [ex_inbound_wellstend; inbound_wellstend];
            
            cntsess=cntsess+1;
            ex_ntrajs_persess(cntsess) = length(inbound_logic);
            n_temp = n_temp + length(inbound_logic);
            
        end
        cnt=cnt+1;
        ex_ntrajs_perday(cnt) = n_temp;
        
    end
    
end

% ---------------------------------------
% Get Moving Average

st=1;
currmax=win;

logic = all_inbound_logic;

% Mean perf vs day
stidx=1;
for i=1:length(ntrajs_perday)
    i;
    perfvsday(i) = mean(logic(stidx:stidx+ntrajs_perday(i)-1));
    stidx = stidx+ntrajs_perday(i);    
end


while currmax<length(logic)
    
    if st<=ceil(win/2),
        perf(st)=mean(logic(1:win));
        currmax=win;
    else
        perf(st)=mean(logic(st-(floor(win/2)):st+(floor(win/2))));
        currmax = st+floor(win/2);
    end
    
    st=st+1;
    
end

%smooth_perf=smooth(perf); % Default is a 5-point moving average
smooth_perf = perf';


% ---------------------------------------
% Get Figure

if figopt==1
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1; 
    end
    
    %Get time (day) axis
    t = 1:length(perf);
    %% Plot data %%
     plot(t, perf',[clr,'.-'], 'Linewidth', 1.1, 'MarkerSize', 16);
    %plot(t, smooth_perf,'k-', 'Linewidth', 3);
    
    %Plot trials in a day
    shift = [0, ntrajs_perday(1:end-1)];
    cumshift = cumsum(shift);
    ntraj_axis_day = ntrajs_perday+cumshift;
    
    xax = [0, ntraj_axis_day];
    for i=1:length(days)
        %line([ntrajs_perday(i) ntrajs_perday(i)], [0 1],'g--');
        plot((ntraj_axis_day(i))*ones(size([0:0.05:1.05])) , [0:0.05:1.05], 'k--');
        currxax = xax(i):xax(i+1);
        % Get triangle on top of shading area
        midpt = floor(length(currxax)/2);
        len1 = midpt;
        spacing1  = (1.05-1)./(len1-1);
        yup1 = 1:spacing1:1.05; % Corrsponds to currax(1):midpt
        len2 = length(currxax) - len1; 
        spacing2 = (1-1.05)./(len2-1);
        yup2 = 1.05:spacing2:1; % Corrsponds to currax(midpt+1):end
        yup = [yup1,yup2];
        
        if mod(i,2)==1 % odd days
            jbfill(currxax,yup, 0*ones(size(currxax)),'k','k',1,1);
        else % even days
            jbfill(currxax,yup, 0*ones(size(currxax)),'k','k',1,1);
        end
    end   
   
    %Backgnd Prob
    %line([1 t(end)], [background_prob  background_prob ], 'k');
    plot(t, 0.5*ones(size(t)),'k', 'Linewidth', 2);
    
    %Title
    title([prefix ' -' 'Inbound Trials (Moving Average-' num2str(win) 'points) Performance'],'FontSize',titlefont,'Fontweight','normal');
    
    % Axis
    xlabmax = xax(end) - mod (xax(end),50); % Get max multiple of 50 that fits into axis
    set(gca,'XTick',[0:50:xlabmax],'XTickLabel',num2str([0:50:xlabmax]'));
    axis([1 t(end)  0 1.1]);
    set(gca,'TickDir','out');
    
    %Axes Names
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal')
    ylabel('Probability of a Correct Response','FontSize',yfont,'Fontweight','normal')
end

if savefig == 1    
    [yyyy,mm,dd] = datevec(date); % if you want to save date in figure name, use this
    figfile = [figdir,prefix,'_in_raw'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end


% % Old
% xax = [cumshift, cumshift(end) + ntrajs_perday(end)];
% for i=1:length(ntrajs_perday)
%     %line([ntrajs_perday(i) ntrajs_perday(i)], [0 1],'g--');
%     plot((ntrajs_perday(i)+cumshift(i))*ones(size([0:0.1:1.1])) , [0:0.1:1.1], 'k--');
%     currxax = xax(i):xax(i+1);
%     if mod(i,2)==1 % odd days
%         jbfill(currxax,1*ones(size(currxax)), 0*ones(size(currxax)),'k','k',1,0.15);
%     else % even days
%         jbfill(currxax,1*ones(size(currxax)), 0*ones(size(currxax)),'k','k',1,0.1);
%     end
%     
% end


