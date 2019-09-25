%    [h] = plotratemap(ratemap, options)
%           options are
%            'fontsize', n    n is the fontsize is point
%            'peakrate', r    set the peak of the colormap to r
%            'minrate', r    do not display ratemaps with peak < r (default 0)
%	     'showmax', 0 or 1  1 indicates that only the maximum rate should
%	     			be labeled
% function [h ch] = plotsj_ratemap(rmap,  varargin)
% 
% minrate = 0;
% peakrate = [];
% fontsize = 12;
% showmax = 0;
% if (~isempty(varargin))
%     assign(varargin{:});
% end
% 
% if (isempty(rmap))
%     return;
% end

function [ch] = plotsj_ratemap_mulepochs(rmap, prefix, day, eps, tet, cell, saveg, varargin)
% Like plotsj_ratemap, but for multiple (right now - 2) epochs
% plotsj_2dmaps also plots separated by trajectory
%  plotsj_ratemap_mulepochs([], 'REf', 8, [2 4], 9, 1, 0);

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
% If rmap is empty and file is specified, go to directory and get data

if (isempty(rmap))
    
    switch prefix
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
    fname = [prefix,'mapfields',daystring];
    load(fname);
    rmap1 = mapfields{day}{eps(1)}{tet}{cell}.smoothedspikerate;
    rmap2 = mapfields{day}{eps(2)}{tet}{cell}.smoothedspikerate;
    % Get everything so you can also look at xticks and yticks if needed
    map1 = mapfields{day}{eps(1)}{tet}{cell};
    map2 = mapfields{day}{eps(2)}{tet}{cell};
    
    %return;
end


% ------------------------------
% Figure and Font Sizes

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/PlaceFields/';
forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

%set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',12);
    tfont = 12; % title font
    xfont = 12;
    yfont = 12;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end

clr = {'b','r','g','c','m','y','k','r'};

% ---------------------------------------


% create a single behav
%cmap = jet(1024) ./ 1.5;
%cmap = hot(1024) ./ 1.5;
%cmap = cmap(100:920,:);

figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end

cmap = jet(1024);
cmap(1,:) = 1;
colormap(cmap);

%Make the subplots for each trajectory and epoch, Normalize color for ratemaps across
% the two epochs, and make sure that the x/y axis is the same for all
% subplots
max1 = max(rmap1(:));
max2 = max(rmap2(:));

% ---------------------------------------
% Set up the bounds to make it look good
bounds = [0 0];
if (isempty(peakrate))
    peakrate = max([max1 max2]);
    bounds(2) = peakrate * 0.65; %0.65;
else
    bounds(2) = peakrate;
end

% Option1: Just set bounds1 to -1
%bounds = [-1 peakrate*.65]; 

% Option2: Tailor bounds(1) to peakrate to get better coloring
if (peakrate > 20)
    minbound = -1;
elseif (peakrate > 10)
    minbound = -.5;
    rmap1(find(rmap1 == -1)) = -.5;
    rmap2(find(rmap2 == -1)) = -.5;
elseif (peakrate > 3)
    minbound = -.1;
    rmap1(find(rmap1 == -1)) = -.1;
    rmap2(find(rmap2 == -1)) = -.1;
else
    minbound = -.01;
    rmap1(find(rmap1 == -1)) = -.01;
    rmap2(find(rmap2 == -1)) = -.01;
end

bounds(1) = minbound;
if (peakrate < .1)
    bounds(2) = 1;
end

% Could also match xticks and yticks to get xbounds and ybounds. You are
% going to get rid of axis labels anyway
x1 = min(map1.xticks);
x2 = min(map2.xticks);
x1m = max(map1.xticks);
x2m = max(map2.xticks);

y1 = min(map1.yticks);
y2 = min(map2.yticks);
y1m = max(map1.yticks);
y2m = max(map2.yticks);

xmin = floor(min([x1 x2]));
xmax = ceil(max([x1m x2m]));
ymin = floor(min([y1 y2]));
ymax = ceil(max([y1m y2m]));


% set up the colormap to be white at some negative values
if (peakrate >= minrate)
    subplot(1,2,1);
    h1 = imagesc(map1.xticks,map1.yticks,rmap1, bounds);
    %h1 = imagesc(flipud(rmap1), bounds); % Produces the same orientation
    % of track as above
    axis([xmin xmax ymin ymax]);
    axis off;
    axis equal;
    title(['Peakrate: ',num2str(floor(max1 * 1)*1)]);
    
    subplot(1,2,2);
    h2 = imagesc(map2.xticks,map2.yticks,rmap2, bounds);
    axis([xmin xmax ymin ymax]);
    axis off;
    axis equal;
    title(['Peakrate: ',num2str(floor(max2 * 1)*1)]);
        
%     ch = colorbar('EastOutside');
%     set(ch, 'FontSize', fontsize);
%     if (showmax)
%         set(ch, 'YTick', floor(peakrate * 0.65), 'YTickLabel', num2str(floor(peakrate * 1))); 
%     end
end


if saveg==1
    figfile = [figdir,prefix,'ratemaps_','d',daystring,'t',tetstring,'c',num2str(cell)];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

