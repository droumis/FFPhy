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

function [h ch] = plotsj_ratemap(rmap, prefix, day, ep, tet, cell, saveg, varargin)
%
%  plotsj_ratemap([], 'REf', 8, 2, 9, 1, 0);,

% This has the rate for all 4 trajectories.

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
    fname = [prefix,'mapfields',daystring];
    load(fname);
    rmap = mapfields{day}{ep}{tet}{cell}.smoothedspikerate;
    %return;
end


% ------------------------------
% Figure and Font Sizes

%figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/PlaceFields/';
figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/Corrln/Egs/';
forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

summdir = figdir;
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
% set up the bounds to make it look good
bounds = [0 0];
if (isempty(peakrate))
    peakrate = max(rmap(:));
    bounds(2) = peakrate * 0.65; %0.65;
else
    bounds(2) = peakrate;
end
if (peakrate > 20)
    minbound = -1;
elseif (peakrate > 10)
    minbound = -.5;
    rmap(find(rmap == -1)) = -.5;
elseif (peakrate > 3)
    minbound = -.1;
    rmap(find(rmap == -1)) = -.1;
else
    minbound = -.01;
    rmap(find(rmap == -1)) = -.01;
end

bounds(1) = minbound;
if (peakrate < .1)
    bounds(2) = 1;
end
% set up the colormap to be white at some negative values
if (peakrate >= minrate)
    h = imagesc(flipud(rmap), bounds);
    ch = colorbar;
    set(ch, 'FontSize', fontsize);
    if (showmax)
        set(ch, 'YTick', floor(peakrate * 0.65), 'YTickLabel', num2str(floor(peakrate * 1))); 
    end
    %title(['Peakrate: ',num2str(floor(peakrate * 1))]);
end

axis off;
axis equal

if saveg==1
    figfile = [figdir,prefix,'ratemap_','d',daystring,'ep',num2str(ep),'t',tetstring,'c',num2str(cell)];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

