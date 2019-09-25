% test to remove artifacts from eeg data
npoints = 15;
thresh = 750;
load /data/loren/shantanu/S1a/EEG/s1aeeg10-4-04.mat
e = eeg{10}{4}{4};
t = geteegtimes(e);
ipoint = zeros(size(e.data));

%load /data/loren/shantanu/S1a/EEG/s1aeeg10-2-04.mat
%e4 = eeg{7}{2}{4};
%high = (abs(e.data) > thresh);
%dh = diff(high);
%d = find(dh == 1);
%for i = 1:length(d)
%    ipoint(round(d(i)-npoints/2):round(d(i)+npoints/2)) = 1;
%end

load /data/loren/shantanu/S1a/s1aDIO10;
pt = DIO{10}{4}{48}.pulsetimes ./ 10000;
eind = lookup(pt(:,1), t);
for i = 1:length(eind)
    ipoint(eind(i):(eind(i)+npoints-1)) = 1;
end

ipoint = logical(ipoint);

newe = interp1(t(~ipoint), e.data(~ipoint), t(ipoint), 'spline');
%newe4 = interp1(t(~ipoint), e4.data(~ipoint), t(ipoint), 'spline')
etmp = e.data;
etmp(ipoint) = newe;


