% test of bessel filter for ripple extraction

%load /data/mkarlsso/Fra/EEG/fraeeg08-2-01
%e = eeg{8}{2}{1};
load /data/mkarlsso/Fra/EEG/fraeeg08-2-29
e = eeg{8}{2}{29};
etimes = geteegtimes(e);

%load /data/mkarlsso/Fra/EEG/fraripple08-2-01
%r = ripple{8}{2}{1};
load /data/mkarlsso/Fra/EEG/fraripple08-2-29
r = ripple{8}{2}{29};

load /data/mkarlsso/Fra/fraspikes08;
s = spikes{8}{2};
%load /data/mkarlsso/Fra/EEG/fraeeg08-3-01
%e = eeg{8}{3}{1};

% create a single mua spikes structure
pushd /data/mkarlsso/frank/frank08
d = dir();
muatimes = [];
for dnum = 3:length(d)
    if (d(dnum).isdir)
	cd(d(dnum).name)
	tmpd = dir();
	load(tmpd(3).name); 
	tmpt = clustdata.params(:,1) / 10000;
	cind = find((tmpt > etimes(1)) & (tmpt < etimes(end)));
	muatimes = [muatimes ; tmpt(cind,1)]; 
	cd ..
	%keyboard
    end
end
popd

muatimes = sort(muatimes);
t = mutimes(1):.01:muatimes(end);
c = histc(muatimes, t);
bar(t,c);

hold on 
plot(allriptimes, ones(size(allriptimes))*20, 'rx', 'MarkerSize', 20);



%load /data/mkarlsso/Fra/EEG/fraripple08-3-01
%r = ripple{8}{3}{1};

load /data/mkarlsso/Fra/fraripples08
%rip = ripples{8}{2}{1}
rip = ripples{8}{2}{29}

% make a variable with all ripple times
allriptimes = rip.midtime;
for t = 1:size(ripples{8}{2})
    if (~isempty(ripples{8}{2}{t}))
	allriptimes = [allriptimes ; ripples{8}{2}{t}.midtime];
    end
end
bins = etimes(1):0.01:etimes(end);
[a b] = flhist(allriptimes, bins);
rfound = find(a);
allrip = b(rfound);

keyboard
;


%plot(e.data(1:10000))

tmpeeg = e.data(1:100000);
tmprip = r.data(1:100000);

samprate = 1500;
lpcutoff = 275;
lpfreq = 250;
hpfreq = 150;
hpcutoff = 125;

% create a butterworth filter

%freq = [lpfreq hpfreq] ./ (samprate / 2);
%[b a] = butter(4, freq);
%load rtripplefilt

%h = spectrum.welch;
%c1psd = psd(h, c1, 'Fs', samprate);
%plot(c1psd)
%c2psd = psd(h, c2, 'Fs', samprate);
%hold on
%plot(c2psd)

load ../Filter/ripplefilter.mat
a = ripplefilter.tf.den;
b = ripplefilter.tf.num;

newa = fliplr(a(2:end));
newb = fliplr(b);

% note that zero padding is essential
%tmpdata = [zeros(100,1); double(e.data(1:100000))];
tmpdata = [zeros(100,1); double(e.data)];

l = length(b);
y = zeros(size(tmpdata));
for i = l:length(tmpdata)
    y(i) = dot(newb, tmpdata(i-18:i)) - dot(newa, y(i-18:i-1));
end

% remove the zero padding
y = y(101:end);


yt = etimes(1:length(y));
plot(yt, y);

% rectify and integrate
%thresh = [50:10:150];
thresh = [4.5:.5:6];
%thresh = [4];
lastval = zeros(20,1);
%lastval = zeros(10,1);
for t = 1:length(thresh)
    tmpy = abs(y);
    v{t} = zeros(size(y));
    m{t} = zeros(size(y));
    s{t} = zeros(size(y));
    rt{t} = [];
    i = 2;
    while (i <= length(y))
	dm = tmpy(i) - m{t}(i-1);
	m{t}(i) = m{t}(i-1) + dm * 0.001;
	ds = abs(tmpy(i) - m{t}(i)) - s{t}(i-1);
	s{t}(i) = s{t}(i-1) + ds * 0.001;
	posgain = mean(lastval);
	lastval(1:19) = lastval(2:20);
	df = tmpy(i) - v{t}(i-1);
%	threshtmp2 = min(m{t}(i) + 2 * s{t}(i), 45);
%	if (tmpy(i) > threshtmp)
	if (df > 0)
	    gain = 1.2;
	    v{t}(i) = v{t}(i-1) + df * posgain;
	else
	    gain = 0.2;
	    v{t}(i) = v{t}(i-1) + df * gain;
        end
	lastval(20) = gain;
	threshtmp = m{t}(i) + thresh(t) * s{t}(i);
	if ((i > 10000) && (v{t}(i) > threshtmp))
	    rt{t} = [rt{t} yt(i)];
	    tmpm = m{t}(i);
	    tmps = s{t}(i);
	    tmpy(i:i+149) = 0;
	    i = i + 149;
	    m{t}(i) = tmpm;
	    s{t}(i) = tmps;
	end
	i = i + 1;
    end
    threshtmp(t) = m{t}(end) + thresh(t) * s{t}(end);
    sprintf('thresh = %d\n', threshtmp(t))
end

ind = 1;
plot(yt, v{ind});
hold on;
plot(yt, y, 'r');
plot(yt, abs(y), 'r');
plot([yt(1) yt(end)], [thresh(ind) thresh(ind)], 'g');
plot(rt{ind}, 120 * ones(size(rt{ind})), 'r.', 'MarkerSize', 20)
plot(rip.midtime, 126 * ones(size(rip.midtime)), 'k.', 'MarkerSize', 20)
plot(allrip, 130 * ones(size(allrip)), 'b.', 'MarkerSize', 20)

plot(yt, r.data(:,3), 'k');
s2 = rip.baseline + 2 * rip.std;
plot([yt(1) yt(end)], [s2 s2], 'k');
