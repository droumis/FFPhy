function [etmp] = sj_removeartifacts1(e, DIO, npoints,figopt)

% test to remove artifacts from eeg data

if nargin<4,
    figopt=0;
end

%npoints = :
%15 for sjstimc
%30 for RE1-Day 1 40 uA

if nargin < 3,
    npoints = 30;   %%% npoints = 15 for SJstimC
    figopt=0;
end

npoints;
thresh = 750;
%load /data25/sjadhav/SJStimC_direct/EEG/sjceeg06-3-02.mat
%e = eeg{6}{3}{2};
t = geteegtimes(e);
ipoint = zeros(size(e.data));

% %load /data/loren/shantanu/S1a/EEG/s1aeeg10-2-04.mat
% %e4 = eeg{7}{2}{4};
% %high = (abs(e.data) > thresh);
% %dh = diff(high);
% %d = find(dh == 1);
% %for i = 1:length(d)
% %    ipoint(round(d(i)-npoints/2):round(d(i)+npoints/2)) = 1;
% %end



%load /data25/sjadhav/SJStimC_direct/sjcDIO06.mat
pt = DIO.pulsetimes ./ 10000;
eind = lookup(pt(1:end,1), t);   

%for i = 1:length(eind)-21   %%% For RE1 Day3 Epoch 2
for i = 1:length(eind)          %%% Can Skip stimulations if we want    
    ipoint(eind(i):(eind(i)+npoints-1)) = 1;
end

ipoint = logical(ipoint);

newe = interp1(t(~ipoint), e.data(~ipoint), t(ipoint), 'spline');
%newe4 = interp1(t(~ipoint), e4.data(~ipoint), t(ipoint), 'spline')
etmp = e.data;
etmp(ipoint) = newe;

%figopt=1;

%%%%%%%%% Plot to check %%%%
if figopt==1,
    
    cnt=0;
    for i =5:length(eind)-5;
        cnt=cnt+1;
        e_stim(cnt,:)=e.data(eind(i+1)-0.5*e.samprate:eind(i+1)+1*e.samprate);
        etmp_stim(cnt,:)=etmp(eind(i+1)-0.5*e.samprate:eind(i+1)+1*e.samprate);
    end
    
    
    figure; hold on;
    redimscreen_land;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    subplot(2,1,1); hold on;
    plot(e_stim(15,:),'k.-','Linewidth',1,'Markersize',6);
    plot(etmp_stim(15,:),'r.','Linewidth',1,'Markersize',6);
    title(['Example raw and filtered EEG around stim']);
    xlim([0 1500]);
    
    subplot(2,1,2); hold on;
    plot(mean(e_stim),'k.-','Linewidth',1,'Markersize',6);
    plot(mean(etmp_stim),'r.','Linewidth',1,'Markersize',6);
    title(['Mean raw and filtered EEG around stim']);
    xlim([0 1500]);
    
end



