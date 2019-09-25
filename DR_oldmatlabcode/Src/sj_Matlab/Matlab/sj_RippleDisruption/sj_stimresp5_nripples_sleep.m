
saveg1=1; 

%%% Compare Sleep Epochs - Number of ripples ###

directoryname = '/data25/sjadhav/RippleInterruption/SJStimC_direct';
prefix = 'sjc';
days = [1:7];  allepochs = [1 3];
clr = {'b','g','c','m','y','k','r'};
tet=6;
std = 5;

n_pre = []; n_post = [];
t_pre = []; t_post = [];

for d = 1:length(days)
    
    day = days(d);
    ripfile = sprintf('%s/%sripples%02dstd%02d.mat', directoryname, prefix, day, std);
    load(ripfile);
    
    for ep = 1:length(allepochs)
        
        epoch = allepochs(ep);
        nrip = length(ripples{day}{epoch}{tet}.startind);
        all_trip = ripples{day}{epoch}{tet}.endtime - ripples{day}{epoch}{tet}.starttime;
        if ep==1
            n_pre = [n_pre; nrip];
            t_pre = [t_pre; 1000*mean(all_trip)];
        end
        if ep==2            
            n_post = [n_post; nrip];
            t_post = [t_post; 1000*mean(all_trip)];
        end
        
    end
    
end

npre5 = n_pre(5); n_pre(5) = []; tpre5 = t_pre(5); t_pre(5) = [];
npost5 = n_post(5); n_post(5) = []; tpost5 = t_post(5); t_post(5) = [];


%%%%%% Plot Nripples  %%%%%%%%
figure; hold on;
redimscreen_halfvert(0);
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');

for i=1:length(n_pre),
    plot([1 2], [n_pre(i) n_post(i)], [clr{i} '.-'], 'Linewidth',2, 'Markersize',12);
end
plot([1 2], [npre5 npost5], ['k.--'], 'Linewidth',2, 'Markersize',12);

plot(1,mean(n_pre),'rs', 'Markersize',12,'Linewidth',4);
plot(2,mean(n_post),'rs', 'Markersize',12, 'Linewidth',4);

preerr = [(mean(n_pre)-sem(n_pre)), (mean(n_pre)+sem(n_pre))];
posterr = [(mean(n_post)-sem(n_post)), (mean(n_post)+sem(n_post))];

ypts = preerr(1):1:preerr(2);
plot(1*ones(size(ypts)),ypts,'r-', 'Linewidth',4);
ypts = posterr(1):1:posterr(2);
plot(2*ones(size(ypts)),ypts,'r-', 'Linewidth',4);

title(['No of ripples (SEP sd:' num2str(std) ') on Tet ' num2str(tet) ' in sleep epochs'],...
    'FontSize',24,'Fontweight','bold');
set(gca,'xtick',[1 2],'xticklabel',{'Pre-Sleep';'Post-Sleep'});
axis([0.8 2.2 0 1400]);
%xlabel('Sleep Epoch');
ylabel('No of ripples');

[h,p,ci] = ttest2(n_pre, n_post, 0.05, 'left' );
text( 1.1, 1000,['p =' num2str(round(p*100)/100)],'FontSize', 24, 'FontWeight','bold');

if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['RippleN_sleep_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['RippleN_sleep_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
end




%%%%%% Plot Tripples  %%%%%%%%
figure; hold on;
redimscreen_halfvert(0);
orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');

for i=1:length(t_pre)
    plot([1 2], [t_pre(i) t_post(i)], [clr{i} '.-'], 'Linewidth',2, 'Markersize',12);
end
plot([1 2], [tpre5 tpost5], ['k.--'], 'Linewidth',2, 'Markersize',12);

plot(1,mean(t_pre),'rs', 'Markersize',12, 'Linewidth',4);
plot(2,mean(t_post),'rs', 'Markersize',12, 'Linewidth',4);

preerr = [(mean(t_pre)-sem(t_pre)), (mean(t_pre)+sem(t_pre))];
posterr = [(mean(t_post)-sem(t_post)), (mean(t_post)+sem(t_post))];

ypts = preerr(1):1:preerr(2);
plot(1*ones(size(ypts)),ypts,'r-', 'Linewidth',4);
ypts = posterr(1):1:posterr(2);
plot(2*ones(size(ypts)),ypts,'r-', 'Linewidth',4);


title(['Ripple duration (SEP sd: ' num2str(std) ') on Tet ' num2str(tet) ' in sleep epochs'],...
    'FontSize',24,'Fontweight','bold');
set(gca,'xtick',[1 2],'xticklabel',{'Pre-Sleep';'Post-Sleep'});
axis([0.8 2.2 100 170]);
%xlabel('Sleep Epoch');
ylabel('Ripple Duration (ms)');
[h,p,ci] = ttest2(t_pre, t_post, 0.05, 'left' );
text( 1.1, 160,['p =' num2str(round(p*100)/100)],'FontSize', 24, 'FontWeight','bold');



if saveg1==1,
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['RippleN_sleep_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
    saveas(gcf,['RippleDur_sleep_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
end



%%%%%%%%%%%%%%% %%%%%%%%


