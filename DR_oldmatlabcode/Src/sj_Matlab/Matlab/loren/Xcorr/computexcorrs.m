function computexcorrs(directoryname, fileprefix, day, bin, tmax)
%function computexcorrs(directoryname, fileprefix, day, bin, tmax)
%
%	Reads in the spike files from the specified day,  
%	computes the unnormalized cross correlations and saves the results for
%	each day.
%
%	assumes spike data stored in the spike file in animdirectory
%
%directoryname - example '/data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix	- folder name where the day's data is stored
%
%day		- the day to process
%
%bin		- the binsize in sec (0.002 recommended)
%
%tmax		- the absolute value of the maximum time for the cross
%		  correlation in sec (1 or 2 recommended)


% load the spike file
load(sprintf('%s/%sspikes%02d.mat', directoryname, fileprefix, day));

% go through every pair of cells and compute the 10 ms lag cross correlation
% out to 1000 ms

xcorr = {};
for d = 1:length(spikes)
    for e = 1:length(spikes{d})
        for t1 = 1:length(spikes{d}{e})
	    for c1 = 1:length(spikes{d}{e}{t1})
		try
		    c1times = spikes{d}{e}{t1}{c1}.data(:,1);
	        catch
		    continue
	        end
		for t2 = t1:length(spikes{d}{e})
		    startcell = 1;
		    if (t2 == t1)
			% start with the current cell to include the
			% autocorrelgram
			startcell = c1;
		    end
		    for c2 = startcell:length(spikes{d}{e}{t2})
			try 
			    c2times = spikes{d}{e}{t2}{c2}.data(:,1);
		        catch
			    continue
		        end
			xcorr{d}{e}{t1}{c1}{t2}{c2} = ...
				spikexcorr(c1times, c2times, 0.002, 2);
		    end
		end
	    end
	end
    end
end
%save(sprintf('%s/%sxcorr%02d.mat', directoryname, fileprefix, day), 'xcorr');
