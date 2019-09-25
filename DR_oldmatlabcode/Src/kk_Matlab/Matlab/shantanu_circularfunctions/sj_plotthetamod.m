function out = sj_plotthetamod(sind, tind, excludetimes, spikes, theta, varargin)
% out = plotthetamod(spike_index, theta_index, excludeperiods, spikes, theta, options)
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
s = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data(:,1);

% get the eeg sample times
t = geteegtimes(theta{tind(1)}{tind(2)}{tind(3)});

tph = theta{tind(1)}{tind(2)}{tind(3)}.data(:,2);
sph = tph(lookup(s, t));



if ~isempty(excludetimes)
    totalexclude = sum(excludetimes(:,2) - excludetimes(:,1));
else
    totalexclude = 0;
end

if ~isempty(s)
    if ~isempty(sph)
        goodspikes = ~isExcluded(s, excludetimes);
    else
        goodspikes = 0;
    end
else
    goodspikes = nan;
end
numgoodspikes = sum(goodspikes);

% This is phase angle at which each spike occurs
% Here, phase is -pi to pi, which is equivalent to 0 to 2*pi
sph = double(sph(goodspikes)) / 10000; 

% Stats and Fit
% ---------------
% Rayleigh and Modulation: Originally in lorenlab Functions folder
stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
[m, ph] = modulation(sph);
phdeg = ph*(180/pi);
% Von Mises Distribution - From Circular Stats toolbox 
[thetahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
thetahat_deg = thetahat*(180/pi);
[prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data


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
figure; hold on;
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
figure;
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
nbins = 50; 
bins = -pi:(2*pi/nbins):pi;
count = histc(sph, bins);

% Make Von Mises Fit
alpha = linspace(-pi, pi, 50)';
[pdf] = circ_vmpdf(alpha,thetahat,kappa);
figure; %redimscreen_2horsubplots
%subplot(1,2,1); 
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

%subplot(1,2,2); hold on;
figure; hold on;

totalspks = sum(count); 
countper = (count./totalspks)*100;
out = bar(bins, countper, 'hist');  set(gca, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
set(gca,'XLim',[-pi pi]); ylabel('% of Spikes'); 
set(out,'FaceColor','r'); set(out,'EdgeColor','r'); 
pdf = pdf.*(countper(binnum)/pdf(binnum));
plot(alpha,pdf,'k','LineWidth',3,'Color','k');


i=1;

%pause
