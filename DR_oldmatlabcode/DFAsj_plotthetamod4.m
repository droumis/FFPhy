function out = DFAsj_plotthetamod4(sind, tind, excludetimes, spikes, theta, varargin)
% out = plotthetamod(spike_index, theta_index, excludeperiods, spikes, theta, options)
% Ver4 : Starting 10Feb2014 - Sync codes with everyone

% For filter framework. See also sj_plotthetamod

% plots the theta moduation histogram for the specified cell based on the
% specified index into the theta structure.
% Excluded time periods are not included in calculation.
%
% Options:
%   'nbins', # of bins (default 24)
%   'appendind', 1 or 0 -- set to 1 to append the cell ind to the
%   output [tetrode cell value].  Default 0.
%
% eg. sind = [ 2 2 18 1] or [2 2 15 1].  tind =[2 2 1]
% sj_plotthetamod([2 2 18 1], [2 2 1], [], spikes, theta);


figopt=0;

nbins = 24;
appendind = 0;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendind'
                appendind = varargin{option+1};
            case 'nbins'
                nbins = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

% get the spike times
if ~isempty(spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data)
    s = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data(:,1);
else
    s=[];
end

% get the eeg sample times and phase
try
    t = geteegtimes(theta{tind(1)}{tind(2)}{tind(3)});
catch
    disp('Stopped in DFAsj_plotthetamod');
    [tind(1) tind(2) tind(3)]
    keyboard; 
   
end
tph = theta{tind(1)}{tind(2)}{tind(3)}.data(:,2);


if ~isempty(s)
    goodspikes = ~isExcluded(s, excludetimes);
else
    goodspikes = [];
end

if length(goodspikes)~=0
    sph = tph(lookup(s, t));  
    sph = double(sph(goodspikes)) / 10000;  % If no spikes, this will be empty
else
    sph = [];
end

numgoodspikes = sum(goodspikes);

% This is phase angle at which each spike occurs
% Here, phase is -pi to pi, which is equivalent to 0 to 2*pi


 
% Stats and Fit
% ---------------
if length(sph)>1
    
    % Rayleigh and Modulation: Originally in lorenlab Functions folder
    stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
    [m, ph] = modulation(sph);
    phdeg = ph*(180/pi);
    % Von Mises Distribution - From Circular Stats toolbox
    [thetahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
    thetahat_deg = thetahat*(180/pi);
    
    
    [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
    
    % Make finer polar plot and overlay Von Mises Distribution Fit.
    % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
    % -------------------------------------------------------------
    nbins = 50;
    bins = -pi:(2*pi/nbins):pi;
    count = histc(sph, bins);
    
    % Make Von Mises Fit
    alpha = linspace(-pi, pi, 50)';
    [pdf] = circ_vmpdf(alpha,thetahat,kappa);
    
else
    
    stats=0;
    m=0;
    phdeg=0;
    kappa=0;
    thetahat_deg=0;
    prayl=0;
    zrayl=0;
    alpha=0;
    pdf=0;
  
end


% Output
% ------

out.index = sind;
out.tetindex = tind;
out.sph = sph;
out.Nspikes = numgoodspikes;
% Output stats also, although you will have to redo this after combining epochs
% Rayleigh test
out.stats = stats;
out.modln = m;
out.phdeg = phdeg;
% Von Mises Fit
out.kappa = kappa;
out.thetahat_deg = thetahat_deg;
out.prayl = prayl;
out.zrayl = zrayl;
out.alpha = alpha;
out.vmpdf = pdf;

    
% if ~isempty(sph)
%     out.sph = sph;   
% else
%     out.sph = zeros(0,1);
% end
    



if (figopt)

    
    % Plot
    % ----
    bins = -pi:(2*pi/nbins):pi;
    count = histc(sph, bins);
    figure; redimscreen_2horsubplots
    subplot(1,2,1);
    hold on;
    out = bar(bins, count, 'hist'); set(gca, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
    set(gca,'XLim',[-pi-0.5 pi]); ylabel('NSpikes'); set(out,'FaceColor','r');
    title(sprintf('cell %d %d %d %d, eegtet %d theta, mod %f, peakph %g, pval %f', sind, tind(3), m, phdeg, stats.p));
    
    % Shantanu - add a plot with %tage of spikes
    % -------------------------------------------
    subplot(1,2,2);
    hold on;
    totalspks = sum(count);
    countper = (count./totalspks)*100;
    out = bar(bins, countper, 'hist');  set(gca, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
    set(gca,'XLim',[-pi-0.5 pi]); ylabel('% of Spikes'); set(out,'FaceColor','r');
    %title(sprintf('cell %d %d %d %d, eegtet %d theta, mod %f, peakph %g, pval %f', sind, tind(3), m, phdeg, stats.p));
    
    % Polar plot  - use rose or polar:
    % --------------------------------
    % [t,r] = rose(sph): input is angle in radians for each spike
    % polar(t,r);
    
    % My range is -pi to pi. This will plot from 0 to 2*pi, and it seems to do this automatically.
    % No need to adjust range %sphn = sph + pi;
    
    [t,r] = rose(sph);
    figure;  redimscreen_figforppt1;
    polar(t,r,'r'); hold on;
    % The peak phase angle
    lims = get(gca,'XLim');
    radius = lims(2);
    xx = radius .* cos(ph); yy = radius .* sin(ph);
    line([0 xx], [0 yy],'LineWidth',4,'Color','k'); %plot(xx,yy,'ko','MarkerSize',4);
    title(sprintf('cell %d %d %d %d, eegtet %d theta, mod %f, peakph %g, pval %f', sind, tind(3), m, phdeg, stats.p));
    
    
    % Make finer polar plot and overlay Von Mises Distribution Fit.
    % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
    % -------------------------------------------------------------
    
    figure; redimscreen_2horsubplots
    subplot(1,2,1);
    hold on;
    out = bar(bins, count, 'hist'); set(gca, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
    set(gca,'XLim',[-pi pi]); ylabel('NSpikes');
    set(out,'FaceColor','r'); set(out,'EdgeColor','r');
    %pdf = pdf.*(max(count)/max(pdf));
    % Instead of maximum - match values at a bin, maybe close to peak
    binnum = lookup(thetahat,alpha);
    pdf = pdf.*(count(binnum)/pdf(binnum));
    plot(alpha,pdf,'k','LineWidth',3,'Color','k');
    title(sprintf('cell %d %d %d %d, eegtet %d theta, Kappa %f, prefang %g, pval %f', sind, tind(3), kappa, thetahat_deg, prayl));
    
    subplot(1,2,2); hold on;
    
    totalspks = sum(count);
    countper = (count./totalspks)*100;
    out = bar(bins, countper, 'hist');  set(gca, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
    set(gca,'XLim',[-pi pi]); ylabel('% of Spikes');
    set(out,'FaceColor','r'); set(out,'EdgeColor','r');
    pdf = pdf.*(countper(binnum)/pdf(binnum));
    plot(alpha,pdf,'k','LineWidth',3,'Color','k');
    title(sprintf('Nspikes %d', totalspks));
    
end


%pause
