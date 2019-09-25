function out = DFAsj_peakdistance_linpos(ind, excludetimes, spikes, linpos, varargin)
% Shantanu - Renamed from calcpairxcorr
% Changing to return only correlation
% peak distance returned by peakdistance_traj or peakdistance_linpos

%function out = calcxcorrmeasures(index, excludetimes, spikes, varargin)
% Calculates the excess correlation and RMS time lag for the specified cell
% pairs using only spikes not excluded by the excludetimes
%
% Shantanu:
% Changing from working for multicellanal to singlecellanal
% Relative Spike Timing vs Peak Distance
% This is also ideal for sequence compression index. See Sens paper - Fig 7
% Need to limit C to 100ms and add a condition of 30 spikes
% Also, could use CoM rather than peak place fields. There must be another function
% calcxcorrmeasure is a better function with more control which also uses
% linfields (trajdata) to calculate overlap rather than getting place fields from scratch
% Only returns time lag - not the excess correlation
% Does not call xcorrms for RMS time lag, rather it is time of max cross-corr
% Place field comparison is distance between peaks
%
% Options:
%   'bin', n 	binsize in sec. Default 0.002 (2 ms)
%   'tmax', n	maximum time for cross correlation. Default 1 sec.


plotfields = 0;
binsize = 2;

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})           
            case 'plotfields'
                plotfields = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

warning('OFF','MATLAB:divideByZero');


tet1=ind(3); cell1=ind(4); tet2=ind(5); cell2=ind(6);

out.index = ind;
out.peakdist = NaN;
out.peakrate1 = NaN;
out.peakrate2 = NaN;


% Check if linpos exists for epoch and then go ahead

if ( (length(linpos{ind(1)}) >= ind(2)) && (isfield(linpos{ind(1)}{ind(2)},'statematrix')) )
    statematrix = linpos{ind(1)}{ind(2)}.statematrix;
    intersections = linpos{ind(1)}{ind(2)}.wellSegmentInfo.distanceToIntersection*1000;
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
     
        try
            spikes1 = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
            spikes2 = spikes{ind(1)}{ind(2)}{ind(5)}{ind(6)}.data(:,1);
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
            out.peakdist = peakdistance;
            out.peakrate1 = peak1;
            out.peakrate2 = peak2;
            
        else
            out.peakdist = NaN;
            out.peakrate1 = NaN;
            out.peakrate2 = NaN;
        end
        
    end
end





