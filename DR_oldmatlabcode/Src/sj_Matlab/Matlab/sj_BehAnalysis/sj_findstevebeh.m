
function [outbound_curve, all_outbound_logic, all_outbound_wellstend, ntrajs_perday] = sj_findstevebeh(animidx, doprobfit)

%% Gets Behavior Curves from Steves summary data
% Shantanu Jadhav, 02/04/11

%  sj_findstevebeh(3, 1);
% Indexes 1-4 are controls M06,M24,M25,M26 and indexes 5-10 are lesion
% animals

if nargin<2,
    doprobfit = 0;
end

dir='/data25/sjadhav/RippleInterruption/Ste';
cd(dir);
load behavperform_2
days=1:10;

animprefix=behavperform(animidx).subject;
prefix=animprefix;

% Outbound
all_outbound_logic = behavperform(animidx).outreward;
outtrials = behavperform(animidx).dayouttrials;

if animidx ==1
    outtrials(1,:)=[0 1];
end
for cntsess=1:20
    ntrajs_persess_out(cntsess)=outtrials(cntsess,2)-outtrials(cntsess,1)+1;
end
ntrajs_perday_out=sum(reshape(ntrajs_persess_out,2,10),1);

%Inbound
all_inbound_logic = behavperform(animidx).inreward;
intrials = behavperform(animidx).dayintrials;

for cntsess=1:20
    ntrajs_persess_in(cntsess)=intrials(cntsess,2)-intrials(cntsess,1)+1;
end
ntrajs_perday_in=sum(reshape(ntrajs_persess_in,2,10),1);


%%% Do Algorithm

if doprobfit == 1
    
    % Outbound
    [pc, lt] = getestprobcorrect(all_outbound_logic, 0.5, 0);
    outbound_curve_witherr = pc(2:end,:);
    outbound_curve = pc(2:end,1);
    outbound_lowerr = pc(2:end,2);
    outbound_uperr = pc(2:end,3);
    outbound_learningtrial = lt;
    
    figure; hold on; redimscreen_portrait;
    %redimscreen;
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
    set(0,'defaultaxeslinewidth',2);
    
    %% Get time (day) axis
    t = 1:size(pc,1)-1;
    %%% Plot data %%
    plot(t', pc(2:end,1),'r.-', 'Linewidth', 2);
    taxis = t';
    data=[pc(2:end,2),pc(2:end,3)];
    jbfill(taxis',pc(2:end,3)',pc(2:end,2)','r','r',1,0.1);
    plot(t', pc(2:end,2),'m--', 'Linewidth', 2);
    plot(t', pc(2:end,3),'m--', 'Linewidth', 2);
    
    %% Plot trials in a day
    shift = [0, ntrajs_perday_out(1:end-1)];
    cumshift = cumsum(shift);
    ntraj_axis_day_out = ntrajs_perday_out+cumshift;
    learning_day_out = min(find(ntraj_axis_day_out>lt));
    
    for i=1:length(days)
        %line([ntrajs_perday(i) ntrajs_perday(i)], [0 1],'g--');
        plot((ntraj_axis_day_out(i))*ones(size([0:0.1:1])) , [0:0.1:1], 'k--');
    end
    
    % Trials in a session
    shifts = [0, ntrajs_persess_out(1:end-1)];
    cumshifts = cumsum(shifts);
    ntraj_axis_sess_out = ntrajs_persess_out+cumshifts;
    
    learning_sess_out = min(find(ntraj_axis_sess_out>lt));
    %% Axis
    axis([1 t(end)  0 1]);
    %axis([1 300  0 1]);
    
    %% Backgnd Prob
    %line([1 t(end)], [background_prob  background_prob ], 'k');
    plot(t, 0.5*ones(size(t)),'k', 'Linewidth', 2);
    
    %% Title
    title([prefix ' - Outbound Trials'],'FontSize',24,'Fontweight','bold');
    
    %% Axes Names
    xlabel('Trial Number')
    ylabel('Probability of a Correct Response')
   
    
    %% Get mean performance as a function of day
    ntraj_start_day_out=[1 ntraj_axis_day_out];
    for i=1:length(days),
        outbound_perf_day(i) = mean(outbound_curve(ntraj_start_day_out(i):ntraj_start_day_out(i+1)));
    end
    
    % if isnan(outbound_learningtrial)
    %     outbound_learningtrial=100;
    %     learning_day=5;
    %     learning_sess=10;
    % end
    
    
    % Inbound
    
    [pc, lt] = getestprobcorrect(all_inbound_logic, 0.5, 0);
    inbound_curve_witherr = pc(2:end,:);
    inbound_curve = pc(2:end,1);
    inbound_lowerr = pc(2:end,2);
    inbound_uperr = pc(2:end,3);
    inbound_learningtrial = lt;
    
    figure; hold on; redimscreen_portrait;
    %redimscreen;
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
    set(0,'defaultaxeslinewidth',2);
    
    %% Get time (day) axis
    t = 1:size(pc,1)-1;
    %%% Plot data %%
    plot(t', pc(2:end,1),'r.-', 'Linewidth', 2);
    taxis = t';
    data=[pc(2:end,2),pc(2:end,3)];
    jbfill(taxis',pc(2:end,3)',pc(2:end,2)','r','r',1,0.1);
    plot(t', pc(2:end,2),'m--', 'Linewidth', 2);
    plot(t', pc(2:end,3),'m--', 'Linewidth', 2);
    
    %% Plot trials in a day
    shift = [0, ntrajs_perday_in(1:end-1)];
    cumshift = cumsum(shift);
    ntraj_axis_day_in = ntrajs_perday_in+cumshift;
    learning_day_in = min(find(ntraj_axis_day_in>lt));
    
    for i=1:length(days)
        %line([ntrajs_perday(i) ntrajs_perday(i)], [0 1],'g--');
        plot((ntraj_axis_day_in(i))*ones(size([0:0.1:1])) , [0:0.1:1], 'k--');
    end
    
    % Trials in a session
    shifts = [0, ntrajs_persess_in(1:end-1)];
    cumshifts = cumsum(shifts);
    ntraj_axis_sess_in = ntrajs_persess_in+cumshifts;
    
    learning_sess_in = min(find(ntraj_axis_sess_in>lt));
    %% Axis
    axis([1 t(end)  0 1]);
    %axis([1 300  0 1]);
    
    %% Backgnd Prob
    %line([1 t(end)], [background_prob  background_prob ], 'k');
    plot(t, 0.5*ones(size(t)),'k', 'Linewidth', 2);
    
    %% Title
    title([prefix ' - Inbound Trials'],'FontSize',24,'Fontweight','bold');
    
    %% Axes Names
    xlabel('Trial Number');
    ylabel('Probability of a Correct Response');
    
     %% Get mean performance as a function of day
    ntraj_start_day_in=[1 ntraj_axis_day_in];
    for i=1:length(days),
        inbound_perf_day(i) = mean(inbound_curve(ntraj_start_day_in(i):ntraj_start_day_in(i+1)));
    end
    
    
    %%%% SAVE
    
    savefile = [animprefix '_outbound'];
    save(savefile,'all_outbound_logic', 'ntrajs_perday_out','outbound_curve_witherr',...
        'outbound_curve','outbound_lowerr','outbound_uperr','outbound_learningtrial',...
        'ntrajs_persess_out','ntraj_axis_day_out','ntraj_axis_sess_out','learning_day_out','learning_sess_out','outbound_perf_day');
    
    savefile = [animprefix '_inbound'];
    save(savefile,'all_inbound_logic', 'ntrajs_perday_in','inbound_curve_witherr',...
        'inbound_curve','inbound_lowerr','inbound_uperr','inbound_learningtrial',...
        'ntrajs_persess_in','ntraj_axis_day_in','ntraj_axis_sess_in','learning_day_in','learning_sess_in','inbound_perf_day');
    
end


