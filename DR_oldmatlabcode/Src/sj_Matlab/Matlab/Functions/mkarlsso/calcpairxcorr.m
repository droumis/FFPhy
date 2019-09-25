function out = calcpairxcorr(ind, excludetimes, spikes, linpos, varargin)
%function out = calcxcorrmeasures(index, excludetimes, spikes, varargin)
% Calculates the excess correlation and RMS time lag for the specified cell
% pairs using only spikes not excluded by the excludetimes
%
%
% Options:
%   'bin', n 	binsize in sec. Default 0.002 (2 ms)
%   'tmax', n	maximum time for cross correlation. Default 1 sec.

bin = 0.0001;
tmax = .5;


for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            
            case 'bin'
                bin = varargin{option+1};
            case 'tmax'
                tmax = varargin{option+1};        
            case 'plotxcorr'
                plotxcorr = varargin{option+1};         
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

if (nargin < 5)
    binsize = 2;
end

warning('OFF','MATLAB:divideByZero');

for i = 1:size(ind,1)
    try
        spikes1 = spikes{ind(1,1)}{ind(1,2)}{ind(i,3)}{ind(i,4)}.data(:,1);
        spikes1 = spikes1(find(~isExcluded(spikes1,excludetimes)));
        spikes2 = spikes{ind(2,1)}{ind(2,2)}{ind(i,5)}{ind(i,6)}.data(:,1);
        spikes2 = spikes2(find(~isExcluded(spikes2,excludetimes)));
    catch
        spikes1 = [];
        spikes2 = [];
    end
    if (~isempty(spikes1) && ~isempty(spikes2))


        %xc = spikexcorr(t1inc, t2inc, bin, tmax);

        tmpcorr = spikexcorr(spikes1, spikes2, bin, tmax);
        maxval = max(tmpcorr.c1vsc2);
        out(i).time = [];
        for addup = 1:maxval
            out(i).time = [out(i).time tmpcorr.time(find(tmpcorr.c1vsc2 == addup))];
        end    
        %out(i).time = tmpcorr.time;
        %out(i).c1vsc2 = tmpcorr.c1vsc2;
        out(i).nspikes1 = length(spikes1);
        out(i).nspikes2 = length(spikes2);
        out(i).ind = ind(i,:);
        out(i).peakdist = [];
        out(i).peakrate1 = [];
        out(i).peakrate2 = [];
    else
        out(i).time = [];
        %out(i).c1vsc2 = [];
        out(i).nspikes1 = [];
        out(i).nspikes2 = [];
        out(i).ind = ind(i,:);
        out(i).peakdist = [];
        out(i).peakrate1 = [];
        out(i).peakrate2 = [];
    end
end

if ( (length(linpos{ind(1,1)}) >= ind(1,2)) && (isfield(linpos{ind(1,1)}{ind(1,2)},'statematrix')) )
    statematrix = linpos{ind(1,1)}{ind(1,2)}.statematrix;
    intersections = linpos{ind(1,1)}{ind(1,2)}.wellSegmentInfo.distanceToIntersection*1000;
    timestep = statematrix.time(2,1) - statematrix.time(1,1);

    %bin the distance from two of the wells, where 0 is the W intersection
    goodtimes = find(~isExcluded(statematrix.time, excludetimes));
    %goodtimes = find(abs(statematrix.linearVelocity(:,2)) > 3);
    behavetimes = statematrix.time(goodtimes);
    
    if ~isempty(goodtimes)
        welldist = statematrix.linearDistanceToWells(goodtimes,1:2);
        bins = [];
        bins{1} = [fliplr([intersections(1):-binsize:0]) [intersections(1)+binsize:binsize:max(welldist(:,1))]]';
        bins{2} = [fliplr([intersections(2):-binsize:0]) [intersections(2)+binsize:binsize:max(welldist(:,2))]]';
        welldistbins = [];
        for i = 1:2
            welldistbins = [welldistbins bins{i}(lookup(welldist(:,i),bins{i}))];
        end

        %get all the unique bins
        [uniquebins, crap, uniquebinindex] = unique(welldistbins,'rows');

        %calculate the distance between each pair of bins
        distancematrix = zeros(size(uniquebins,1),size(uniquebins,1));
        for i = 1:size(uniquebins,1)
            tmp = repmat(uniquebins(i,:),size(uniquebins,1),1);
            distancematrix(:,i) = max(abs(tmp-uniquebins),[],2);
        end

        %create a smoothing matrix based on a gaussian curve
        distancematrix2 = round(distancematrix/binsize)+1;
        smoothcurve = gaussian(2,max(distancematrix2(:))*2);
        smoothcurve = smoothcurve(round(length(smoothcurve)/2):end);
        smoothcurve(end+1) = 0;
        distancematrix2 = smoothcurve(distancematrix2);
        distancematrix2 = distancematrix2 * diag(1./sum(distancematrix2)); %normalize each bins smoothing to 1

        %calculate the smoothed occupancy
        occupancy = rowcount(uniquebins,welldistbins)*timestep;
        smoothoccupancy = distancematrix2'*occupancy;

        for i = 1:size(ind,1)

            try
                spikes1 = spikes{ind(1,1)}{ind(1,2)}{ind(i,3)}{ind(i,4)}.data(:,1);
                spikes2 = spikes{ind(2,1)}{ind(2,2)}{ind(i,5)}{ind(i,6)}.data(:,1);
            catch
                spikes1 = [];
                spikes2 = [];
            end
            if (~isempty(spikes1) && ~isempty(spikes2))
                %keep the time and the posindex for each valid spike
                spikes1 = spikes1(~isExcluded(spikes1,excludetimes));
                spikes2 = spikes2(~isExcluded(spikes2,excludetimes));

                %calculate the smoothed spike counts
                spikecount1 = rowcount(uniquebins,uniquebins(uniquebinindex(lookup(spikes1,behavetimes)),:));
                smoothspikecount1 = distancematrix2'*spikecount1;
                spikecount2 = rowcount(uniquebins,uniquebins(uniquebinindex(lookup(spikes2,behavetimes)),:));
                smoothspikecount2 = distancematrix2'*spikecount2;

                %calculate occupancy normalized firing rates
                occnormrate1 = smoothspikecount1./smoothoccupancy;
                occnormrate2 = smoothspikecount2./smoothoccupancy;

                %calculate the distance bateen the peaks
                [peak1, peakind1] = max(occnormrate1);
                [peak2, peakind2] = max(occnormrate2);
                peakdistance = distancematrix(peakind2,peakind1);
                %peakdistance = normoverlap(occnormrate1,occnormrate2);
                out(i).peakdist = peakdistance;
                out(i).peakrate1 = peak1;
                out(i).peakrate2 = peak2;

            else
                out(i).peakdist = [];
                out(i).peakrate1 = [];
                out(i).peakrate2 = [];
            end
        end
    end
end


        


