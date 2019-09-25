
% Similar to plotsj_ratemap
% Load linfields and plot linearized rate for inbond and outbound trajctories separately


function [h ch] = plotsj_linearrate(trajdata, prefix, day, ep, tet, cell, saveg, varargin)
%
%  plotsj_linearrate([], 'HPa', 8, 2, 17, 5, 0);


minrate = 0;
peakrate = [];
fontsize = 12;
showmax = 1;
if nargin<2,
    prefix=[];
end
if nargin<3,
    day=[];
end
if nargin<4,
    ep=[];
end
if nargin<5,
    tet=[];
end
if nargin<6,
    cell=[];
end
if nargin<7,
    saveg=[];
end
if (~isempty(varargin))
    assign(varargin{:});
end

%-------------------------------------------------------
% If trajdata is empty and file is specified, go to directory and get data

if (isempty(trajdata))
    
    switch prefix
        case 'HPa'
            directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
        case 'HPb'
            directoryname = '/data25/sjadhav/HPExpt/HPb_direct/';
        case 'HPc'
            directoryname = '/data25/sjadhav/HPExpt/HPc_direct/';
        case 'RE1'
            directoryname = '/data25/sjadhav/RippleInterruption/RE1_direct';
        case 'RCa'
            directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct';
        case 'RCb'
            directoryname = '/data25/sjadhav/RippleInterruption/RCb_direct';
        case 'RCc'
            directoryname = '/data25/sjadhav/RippleInterruption/RCc_direct';
        case 'RCd'
            directoryname = '/data25/sjadhav/RippleInterruption/RCd_direct';
        case 'REc'
            directoryname = '/data25/sjadhav/RippleInterruption/REc_direct';
        case 'REd'
            directoryname = '/data25/sjadhav/RippleInterruption/REd_direct';
        case 'REe'
            directoryname = '/data25/sjadhav/RippleInterruption/REe_direct';
        case 'REf'
            directoryname = '/data25/sjadhav/RippleInterruption/REf_direct';
    end
    
    cd(directoryname);

    if (day < 10)
        daystring = ['0',num2str(day)];
    else
        daystring = num2str(day);
    end
    if (tet < 10)
        tetstring = ['0',num2str(tet)];
    else
        tetstring = num2str(tet);
    end
    fname = [prefix,'linfields',daystring];
    load(fname);
    trajdata = linfields{day}{ep}{tet}{cell};
    %return;
end


% ------------------------------
% Figure and Font Sizes

figdir = '/data25/sjadhav/HpExpt/Figures/PlaceFields/';
forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

summdir = figdir;
%set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',30);
tfont = 40;
xfont = 30;
yfont = 30;

clr = {'k','r','g','c','m','y','b','r'};

% ---------------------------------------



% Outbound
% --------
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
idx=[1 3];
binsize = 2; %cm
for i=1:length(idx)
    curridx = idx(i);
    xaxis = binsize*(1:size(trajdata{curridx},1));
    plot(xaxis,trajdata{curridx}(:,5),[clr{i} '-'],'LineWidth',4);  
end
title('Outbound Trajectories','FontSize',40,'Fontweight','normal');
xlabel ('Position along linear trajectory (cm)','FontSize',30,'Fontweight','normal');
ylabel ('Firing Rate (Hz)','FontSize',30,'Fontweight','normal');
legend('Outbound Left','Outbound Right');

% Inbound
% --------
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
idx=[2 4];
binsize = 2; %cm
for i=1:length(idx)
    curridx = idx(i);
    xaxis = binsize*(1:size(trajdata{curridx},1));
    xaxisr = xaxis(end:-1:1);
    ydata = trajdata{curridx}(:,5);
    ydatar = ydata(end:-1:1);
    plot(xaxis,ydatar,[clr{i} '-'],'LineWidth',4);  
end
set(gca,'XTick',[0:50:200],'XTickLabel',[200:-50:0]');
title('Inbound Trajectories','FontSize',40,'Fontweight','normal');
xlabel ('Position along linear trajectory (cm)','FontSize',30,'Fontweight','normal');
ylabel ('Firing Rate (Hz)','FontSize',30,'Fontweight','normal');
legend('Inbound Left','Inbound Right');



%axis off;
%axis equal

if saveg==1
    figfile = [figdir,prefix,'linearrate_','d',daystring,'ep',num2str(ep),'t',tetstring,'c',num2str(cell)];
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end

