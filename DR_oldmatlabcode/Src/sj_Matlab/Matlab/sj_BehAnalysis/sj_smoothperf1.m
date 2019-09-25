
function [perf, smooth_perf, perfvsday] = sj_smoothperf1 (prefix, win, logic, ntrajs_perday, type, figopt)
% eg [perf, smooth_perf, perfvsday] = sj_smoothperf1 ('REe', 10, all_outbound_logic, ntrajs_perday_out, 'Out',0)

%load REd_outbound
%prefix='REd';
% win=10;
st=1;
currmax=win;


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

if figopt==1
    figure; hold on; redimscreen_figforppt1;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
    set(0,'defaultaxeslinewidth',2);
    
    %Get time (day) axis
    t = 1:length(perf);
    %% Plot data %%
    plot(t, perf','r-', 'Linewidth', 2.5);
    %plot(t, smooth_perf,'k-', 'Linewidth', 3);
    
    %Plot trials in a day
    shift = [0, ntrajs_perday(1:end-1)];
    cumshift = cumsum(shift);
    xax = [cumshift, cumshift(end) + ntrajs_perday(end)];
    for i=1:length(ntrajs_perday)
        %line([ntrajs_perday(i) ntrajs_perday(i)], [0 1],'g--');
        plot((ntrajs_perday(i)+cumshift(i))*ones(size([0:0.1:1.1])) , [0:0.1:1.1], 'k--');
        currxax = xax(i):xax(i+1);
        if mod(i,2)==1 % odd days
            jbfill(currxax,1*ones(size(currxax)), 0*ones(size(currxax)),'k','k',1,0.15);
        else % even days
            jbfill(currxax,1*ones(size(currxax)), 0*ones(size(currxax)),'k','k',1,0.1);
        end
            
    end
    %Axis
    axis([1 t(end)  0 1.05]);
    %Backgnd Prob
    %line([1 t(end)], [background_prob  background_prob ], 'k');
    plot(t, 0.5*ones(size(t)),'k', 'Linewidth', 2);
    %Title
    title([prefix ' -' type 'bound Trials (Moving Average-' num2str(win) 'points) Performance'],'FontSize',20,'Fontweight','bold');
    %Axes Names
    xlabel('Trial Number')
    ylabel('Probability of a Correct Response')
end