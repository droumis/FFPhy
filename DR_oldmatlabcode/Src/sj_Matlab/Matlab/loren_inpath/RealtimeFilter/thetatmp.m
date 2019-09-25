% plot theta from all electrodes 
adir = '/data/mkarlsso/Bon/'
aname = 'bon'
d = 4;
e = 2;

c = loaddatastruct(adir, aname, 'cellinfo');
task = loaddatastruct(adir, aname, 'task');
for t = 1:length(c{d}{e})
    loadstr = sprintf('load %sEEG/%stheta%02d-%d-%02d', adir, aname, d, e, t);
    try
	eval(loadstr);
    catch
	continue;
    end
    ttmp = theta{d}{e}{t}.data;
    if (length(c{d}{e}{t})) 
	if (isfield(c{d}{e}{t}{1}, 'area'))
	    if (strcmp(c{d}{e}{t}{1}.area, 'CA3'))
		plot(ttmp(10000:20000), 'b');
		hold on;
	    elseif (strcmp(c{d}{e}{t}{1}.area, 'CA1'))
		plot(ttmp(10000:20000), 'r');
		hold on;
	    end
	    title(sprintf('exposure %d', task{d}{e}.exposure));
	end
    end
end
