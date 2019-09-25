function sj_getcellpoplninfo(prefix, varargin)
% INCOMPLETE - Switch to using datafilter framework. Easier to combine across epochs when there is a mismatch.
% - This format is good for updating tags. Do Later


% Shantanu - Jan 2013. Use cellinfo structure to getinfo about fir rate,
% spike width, etc for cell population for given days
% eg. sj_getcellpoplninfo('HPa','days', 1:4);

% See also getcellprop.m: it calls "computecsi" and "modulation.m"

% Can be used in conjunction with sj_plotspike_waveform, which can plot
% spike shapes. Also see spikewaveforms.m in kkay/ folder - faster

%---------------------------------------------------------


if nargin<1,
    keyboard
    error('Please enter Expt Prefix at least');
end
% assign the options - Can also use procOptions
days = 1:8;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'days'
            days = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end  
if min(days)==1 && max(days)==8
    disp('Using days 1:8 - default');
end

% Variables
plotfrline1 = 0.1; % Hz Lines to draw along firing rate axes
plotfrline2 = 7;   % Hz


% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/';
%figdir = '/data25/sjadhav/HPExpt/Figures/CellDescription/';
datadir = '/data25/sjadhav/';
summdir = figdir;

set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end
clr = {'b','r','g','c','m','y','k','r'};

% ---------------------------------------


% -----------------------------------------
% SET DATA
% -------------------------------------------

switch prefix
    case 'HPa'
        rawdir = '/data25/sjadhav/HPExpt/HPa/';
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
    case 'HPb'
        rawdir = '/data25/sjadhav/HPExpt/HPb/';
        directoryname = '/data25/sjadhav/HPExpt/HPb_direct/';
end


clr = {'b','g','c','m','y','k','r'};

% --------------------------------------------------
%  %Load the cellinfo and task file
% -------------------------------------------------

cellinfofile = sprintf('%s/%scellinfo.mat', directoryname, prefix);
load(cellinfofile);




% Loop over days and gather cells - firing rates, spike width, propbursts, numspikes, csi(?) and tags if they exist

for d=1:length(days)
    currday = days(d);
    
    % Load task file for day
    taskfile = sprintf('%s/%stask%02d.mat', directoryname, prefix,currday);
    load(taskfile);
    
    for ep = 1:length(cellinfo{currday})
        currtype = task{currday}{ep}.type;
        if ~strcmp(currtype(1:5),'sleep') % only run epochs
            for tet = 1:length(cellinfo{currday}{ep})
                if ~isempty(cellinfo{currday}{ep}{tet})
                    for cell =1:length(cellinfo{currday}{ep}{tet})
                        if ~isempty(cellinfo{currday}{ep}{tet})
                            
                            fr(currday,ep,tet,cell)=cellinfo{currday}{ep}{tet}{cell}.meanrate;
                            
                        end
                    end
                end
            end
        end
    end                 
                            
    
    
    
end




    

    

    

% 
% %% ------------------------------------------------
% % PLOT
% % -------------------------------------------------
% 
% % All Rewarded vs Not Rewarded
% % ----------------------------
% 
% figure; hold on;
% %redimscreen_figforppt1;
% %set(gcf,'Position',[0 Screen(4)*0.55 Screen(3)*0.2 Screen(4)*0.4])
% %redimscreen_2versubplots
% set(gcf,'Position',[16 44 900 1045]);
% 
% % Raster
% % ---------
% subplot(3,1,[1 2]); hold on; % Raster - rewarded and then unrewarded
% for c=1:length(allrew_spks)
%     tmps = allrew_spks{c};
%     plotraster(tmps,c*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','b');
% end
% push = length(allrew_spks);
% for c=1:length(allnorew_spks)
%     tmps = allnorew_spks{c};
%     plotraster(tmps,(push+c)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','r');
% end
%     
% set(gca,'XLim',[-pret+300 postt-300]);
% set(gca,'YLim',[0 push+length(allnorew_spks)+1]);
% set(gca,'XTick',[]); set(gca,'YTick',[]);
% ypts = 0:1:push+length(allnorew_spks)+1;
% xpts = 0*ones(size(ypts));
% plot(xpts , ypts, 'k--','Linewidth',2);



% ----------------------------------------------------------------------------






i=1;
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



