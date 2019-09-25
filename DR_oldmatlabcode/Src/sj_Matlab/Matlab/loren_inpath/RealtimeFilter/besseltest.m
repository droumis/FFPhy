% test of bessel filter for theta phase extraction

load /data/loren/mkarlsso/Fra/EEG/fraeeg08-4-01
e = eeg{8}{4}{1};

%plot(e.data(1:10000))

tmpdata = e.data(1:10000);

samprate = 1500;
lpfreq = 12;
hpfreq = 6;

12 / 750 * pi;
 

% first create lowpass filter at 12 Hz
[b, a] = besself(4, lpfreq * 2 * pi);
[newblp, newalp] = bilinear(b, a, samprate);

% now create lowpass filter at 6 Hz for subtraction
[b, a] = besself(4, lpfreq * 2 * pi);
[newbhp, newahp] = bilinear(b, a, samprate);

out = filter(newblp, newalp, tmpdata);
out2 = filter(newbhp, newahp, out);

r = out - out2;

%save theta.mat new* out* r*


% apply the low pass filter
tmpdata = e.data(1:100000);
y = zeros(size(tmpdata));
b = fliplr(newblp);
a = fliplr(newalp);
a = a(1:(end-1));
for i = 5:length(tmpdata)
    y(i) = dot(b, tmpdata(i-4:i)) - dot(a, y(i-4:i-1));
end

% apply the high pass filter
tmpy = zeros(size(tmpdata));
b = fliplr(newbhp);
a = fliplr(newahp);
a = a(1:(end-1));
for i = 5:length(tmpdata)
    tmpy(i) = dot(b, y(i-4:i)) - dot(a, tmpy(i-4:i-1));
end

y = y - tmpy;

% 77 degree phase lag from xcorr of y and tmpdata (1.35 rad)

% 210 samples per wave, 45 sample offset



% predict timing of next peak:

y = zeros(size(tmpdata));
tmpy = zeros(size(tmpdata));
s = zeros(size(tmpdata));
by = fliplr(newblp);
ay = fliplr(newalp);
ay = ay(1:(end-1));
btmpy = fliplr(newbhp);
atmpy = fliplr(newahp);
atmpy = atmpy(1:(end-1));
epeakind = [];
peakind = [];
troughind = [];
nepeaks = 1;
npeaks = 1;
ntroughs = 1;
minspacing = 1500 / 12; % 12 Hz spacing
for i = 5:length(tmpdata)
    y(i) = dot(by, tmpdata(i-4:i)) - dot(ay, y(i-4:i-1));
    tmpy(i) = dot(btmpy, y(i-4:i)) - dot(atmpy, tmpy(i-4:i-1));
    s(i) = y(i) - tmpy(i);
    % check for trough
    if ((s(i-1) < s(i-2)) & (s(i) >= s(i-1)))
	if (ntroughs < 3)
	    troughind(ntroughs) = i;
	    ntroughs = ntroughs + 1;
	elseif (i - troughind(ntroughs-1) > minspacing)
	    troughind(ntroughs) = i;
	    ntroughs = ntroughs + 1;
	    % trough found
	    % extrapolate to next peak
	    if (ntroughs > 2) 
		%epeakind(nepeaks) = i + (troughind(ntroughs-1) - ...
		%troughind(ntroughs-2))/2 - 45;
		epeakind(nepeaks) = i + ...
		    mean(diff(troughind(ntroughs-3:ntroughs-1)))/2 - 45;
	    else
		epeakind(nepeaks) = i + 105 - 45;
	    end
	    nepeaks = nepeaks + 1;
	end
    end
end

