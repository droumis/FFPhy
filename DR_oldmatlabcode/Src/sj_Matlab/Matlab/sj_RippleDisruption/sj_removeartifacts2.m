function [etmp] = sj_removeartifacts2(rip, column, DIO, npoints);

% test to remove artifacts from eeg data

if nargin < 3,
    npoints = 30;
end

thresh = 750;
%load /data25/sjadhav/SJStimC_direct/EEG/sjceeg06-3-02.mat
%e = eeg{6}{3}{2};
t = geteegtimes(rip);
ripdata = double(rip.data(:,column)); 
ipoint = zeros(size(ripdata));

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
eind = lookup(pt(:,1), t);
for i = 1:length(eind)
    ipoint(eind(i):(eind(i)+npoints-1)) = 1;
end

ipoint = logical(ipoint);
%ipoint = logical(ipoint(1:length(ripdata)));

newe = interp1(t(~ipoint), ripdata(~ipoint), t(ipoint), 'spline');
%newe4 = interp1(t(~ipoint), e4.data(~ipoint), t(ipoint), 'spline')
etmp = ripdata;
etmp(ipoint) = newe;


