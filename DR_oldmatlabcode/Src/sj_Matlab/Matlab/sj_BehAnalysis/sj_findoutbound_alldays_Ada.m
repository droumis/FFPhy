
function [outbound_curve, all_outbound_logic, all_outbound_wellstend, ntrajs_perday] = sj_findoutbound_alldays_Ada (dir, days,epochs,prefix, doprobfit, days_ex, savedata, savefig)

% 12/27/12
% sj_findoutbound_alldays_Ada ('/data25/sjadhav/HPExpt/Ada_direct',[1],[1:8],'Ada',1);
% For single day expt - divide into sessions

%% Calls sj_day_findoutbound by looping over all days, and returns vector of
% outbound_logic (Correct incorrect), wells exited and entered, and number
% of outbound trajectories per day
% Shantanu Jadhav, 12/30/09
% sj_findoutbound_alldays ('/data25/sjadhav/RippleInterruption/REe_direct',[1:8],[2 4],'REe', 1,[9 10],1,0);

% 02/28/11, Shantanu Jadhav
% Add days_ex: to restart algorithm for extra days at end with new
% conditions
% sj_findoutbound_alldays ([1:8],[2 4],'REd', 1,[9 10]);
% sj_findoutbound_alldays ([1:8],[2 4],'REf', 1,[9 10] );
% sj_findoutbound_alldays ([1:9],[],'Cor', 1);


if nargin<5,
    doprobfit = 0;
end

if nargin<6,
    days_ex=[];
end

if nargin<7,
    savedata=0;
end

if nargin<8,
    savefig=0;
end

% ------------------------------
% Figure and Font Sizes

clr = 'r'; % Main line r/b/k
clr2 = 'r'; % Confidence bounds m/b/k

forppr = 1; 
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Ada_direct/Figures/';
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

epflag=0; % If epoch is not specified, get run epochs from task file
i=1;
%%Initialize
all_outbound_logic = [];
all_outbound_wellstend = [];
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
    
    day;
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
        [outbound_logic, outbound_wellstend] = sj_day_findoutbound (linpos,day,epoch);
        
        %% Add to cumulative count
        all_outbound_logic = [all_outbound_logic; outbound_logic];
        all_outbound_wellstend = [all_outbound_wellstend; outbound_wellstend];
        
        cntsess=cntsess+1;
        ntrajs_persess(cntsess) = length(outbound_logic);
        n_temp = n_temp + length(outbound_logic);
        
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
    ex_outbound_logic = [];
    ex_outbound_wellstend = [];
    
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
            [outbound_logic, outbound_wellstend] = sj_day_findoutbound (linpos,day,epoch);
            
            %% Add to cumulative count
            ex_outbound_logic = [ex_outbound_logic; outbound_logic];
            ex_outbound_wellstend = [ex_outbound_wellstend; outbound_wellstend];
            
            cntsess=cntsess+1;
            ex_ntrajs_persess(cntsess) = length(outbound_logic);
            n_temp = n_temp + length(outbound_logic);
            
        end
        cnt=cnt+1;
        ex_ntrajs_perday(cnt) = n_temp;
        
    end
    
end



%%% Do Algorithm

if doprobfit == 1
    
    [pc, lt] = getestprobcorrect(all_outbound_logic, 0.5, 0);
    outbound_curve_witherr = pc(2:end,:);
    outbound_curve = pc(2:end,1);
    outbound_lowerr = pc(2:end,2);
    outbound_uperr = pc(2:end,3);
    outbound_learningtrial = lt;
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1; 
    end
    %redimscreen;
    %orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
     
    %% Get time (day) axis
    t = 1:size(pc,1)-1;
    %%% Plot data %%
    plot(t', pc(2:end,1),clr, 'Linewidth', 2);
    taxis = t';
    data=[pc(2:end,2),pc(2:end,3)];
    jbfill(taxis',pc(2:end,3)',pc(2:end,2)',clr,clr,1,1);
    plot(t', pc(2:end,2),[clr2,'--'], 'Linewidth', 2);
    plot(t', pc(2:end,3),[clr2,'--'], 'Linewidth', 2);
    
    %% Plot trials in a day
    shift = [0, ntrajs_persess(1:end-1)];
    cumshift = cumsum(shift);
    ntraj_axis_sess = ntrajs_persess+cumshift;
    learning_day = min(find(ntraj_axis_sess>lt));
    
    xax = [0, ntraj_axis_sess];
    for i=1:cntsess
        %line([ntrajs_perday(i) ntrajs_perday(i)], [0 1],'g--');
        plot((ntraj_axis_sess(i))*ones(size([0:0.05:1.05])) , [0:0.05:1.05], 'k--');
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
    
    % Trials in a session
    shifts = [0, ntrajs_persess(1:end-1)];
    cumshifts = cumsum(shifts);
    ntraj_axis_sess = ntrajs_persess+cumshifts;
    
    learning_sess = min(find(ntraj_axis_sess>lt));
       
    %% Backgnd Prob
    %line([1 t(end)], [background_prob  background_prob ], 'k');
    plot(t, 0.5*ones(size(t)),'k', 'Linewidth', 2);
    
    %% Title
    title([prefix ' - Outbound Trials'],'FontSize',titlefont,'Fontweight','normal');
    
    %% Axes Names
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal')
    ylabel('Probability of a Correct Response','FontSize',yfont,'Fontweight','normal')
    
     %% Axis
    xlabmax = xax(end) - mod (xax(end),50); % Get max multiple of 50 that fits into axis
    set(gca,'XTick',[0:50:xlabmax],'XTickLabel',num2str([0:50:xlabmax]'));
    axis([1 t(end)  0 1.1]);
    set(gca,'TickDir','out');
    
    %% Get mean performance as a function of day
    ntraj_start_day=[1 ntraj_axis_sess];
    for i=1:length(days),
        outbound_perf_day(i) = mean(outbound_curve(ntraj_start_day(i):ntraj_start_day(i+1)));
    end
    
    
end

if savefig == 1    
    [yyyy,mm,dd] = datevec(date); % if you want to save date in figure name, use this
    figfile = [figdir,prefix,'_out'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end
    

if isnan(outbound_learningtrial)
    %outbound_learningtrial=ntraj_axis_day(end);
    if isempty(find(outbound_lowerr<0.5, 1))
        outbound_learningtrial=1;
        learning_day=1;
        learning_sess=1;
    else
        outbound_learningtrial=length(all_outbound_logic);
        learning_day=8;
        learning_sess=16;
    end
end

ntrajs_perday_out=ntrajs_perday; ntrajs_persess_out=ntrajs_persess;

%% If you have extra days, make new plot
ex_outbound_perf_day=[];  

if ~isempty(days_ex)
    [pc2, lt2] = getestprobcorrect(ex_outbound_logic, 0.5, 2);
    ex_outbound_curve_witherr = pc(2:end,:);
    ex_outbound_curve=pc2(2:end,1);
    npc=[pc;pc2];
    noutbound_curve = npc(2:end,1);
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1; 
    end
    %redimscreen;
    %orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
   
     %% Get time (day) axis
    t = 1:size(npc,1)-1;
    %%% Plot data %%
    plot(t', npc(2:end,1),clr, 'Linewidth', 2);
    taxis = t';
    jbfill(taxis',npc(2:end,3)',npc(2:end,2)',clr,clr,1,0.1);
    plot(t', npc(2:end,2),[clr2,'--'], 'Linewidth', 2);
    plot(t', npc(2:end,3),[clr2,'--'], 'Linewidth', 2);
    
    n_ntrajs_perday = [ntrajs_perday, ex_ntrajs_perday];
    
    %% Plot trials in a day
    shift = [0, n_ntrajs_perday(1:end-1)];
    cumshift = cumsum(shift);
    n_ntraj_axis_day = n_ntrajs_perday+cumshift;
    
    for i=1:length(days)+length(days_ex)
        %line([ntrajs_perday(i) ntrajs_perday(i)], [0 1],'g--');
        plot((n_ntraj_axis_day(i))*ones(size([0:0.05:1.05])) , [0:0.05:1.05], 'k--');
    end
    
    plot((n_ntraj_axis_day(length(days)))*ones(size([0:0.05:1.05])) , [0:0.05:1.05], 'c--','Linewidth',3);
    
      %% Axis
    axis([1 t(end)  0 1]);
    
    %% Backgnd Prob
    %line([1 t(end)], [background_prob  background_prob ], 'k');
    plot(t, 0.5*ones(size(t)),'k', 'Linewidth', 2);
    
        
    %% Title
    title([prefix ' - Outbound Trials+Extra Days'],'FontSize',titlefont,'Fontweight','normal');
    
    %% Axes Names
    xlabel('Trial Number','FontSize',xfont,'Fontweight','normal')
    ylabel('Probability of a Correct Response','FontSize',yfont,'Fontweight','normal')
    
     %% Get mean performance per day for Extra Days 
    n_ntraj_start_day=[1 n_ntraj_axis_day];
    for i=1:length(days_ex),
        ex_outbound_perf_day(i) = mean(noutbound_curve(n_ntraj_start_day(length(days)+i):n_ntraj_start_day(length(days)+i+1)));
    end
    
else
    
    ex_outbound_logic=[];
    ex_outbound_curve=[];
    ex_outbound_curve_witherr=[];
    
end  %% days_ex

if savedata == 1
    
    savefile = [prefix '_outbound'];
    save(savefile,'all_outbound_logic', 'all_outbound_wellstend', 'ntrajs_perday_out','outbound_curve_witherr',...
        'outbound_curve','outbound_lowerr','outbound_uperr','outbound_learningtrial',...
        'ntrajs_persess_out','ntraj_axis_day','ntraj_axis_sess','learning_day','learning_sess','outbound_perf_day',...
        'ex_outbound_perf_day','ex_ntrajs_persess','ex_ntrajs_perday','ex_outbound_logic','ex_outbound_curve','ex_outbound_curve_witherr');
end

keyboard;
%cd(currdir);




%savefile = [prefix '_outbound'];
%save(savefile,'all_outbound_logic',
%'all_outbound_wellstend','ntrajs_perday','outbound_curve');