function out = DFAsj_peakdistance_traj(index, excludetimes, linfields, varargin)
% Shantanu - Derived from calcpairxcorr - this does only peak distance
% based on trajdata in linfields, not xcorr
%
%function out = calcxcorrmeasures(index, excludetimes, spikes, varargin)
% Calculates the excess correlation and RMS time lag for the specified cell
% pairs using only spikes not excluded by the excludetimes
%
% Shantanu:
% Changing peak distance to be done using trajdata in linfields
% Correlation is returned by DFSsj_calcpairxcorr
%
% Changing from working for multicellanal to singlecellanal
% Relative Spike Timing vs Peak Distance
% This is also ideal for sequence compression index. 
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

binsize = 2; % 2cm for trajdata - default
thresh = 3;

if ~isempty(excludetimes)
    excludetimes = [];
end

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})            
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


tet1=index(3); cell1=index(4); tet2=index(5); cell2=index(6);

out.index = index;
out.peakdist = NaN;
out.peakrate1 = NaN;
out.peakrate2 = NaN;
    
% Peak distance based on linfields
trajdata1 = linfields{index(1)}{index(2)}{index(3)}{index(4)}; %All trajectories
trajdata2 = linfields{index(1)}{index(2)}{index(5)}{index(6)};


%peaks = []; 
peaks1 = []; peaks2 = [];
for traj = 1:length(trajdata1)
    if length(trajdata1)>=traj && length(trajdata2)>=traj
        %if ~isempty(trajdata1{traj}) && ~isempty(trajdata2{traj})
        %peaks = [peaks max([max(trajdata1{traj}(:,5)) max(trajdata2{traj}(:,5))])]; % This is combining both cells
        peaks1 = [peaks1 (max(trajdata1{traj}(:,5)))];
        peaks2 = [peaks2 (max(trajdata2{traj}(:,5)))];
        %end
    end
end


[peakrate1,peaktraj1] = max(peaks1);
[peakrate2,peaktraj2] = max(peaks2);

if peaktraj1(1)==peaktraj2(1),
    maxtraj = peaktraj1(1);
else
    maxtraj = [peaktraj1(1) peaktraj2(1)];
end

 if length(maxtraj)==1  % Same maximum trajectory - Get peak distance
     
     out.peakrate1 = nanmax(trajdata1{maxtraj}(:,5));
     out.peakrate2 = nanmax(trajdata2{maxtraj}(:,5));
     
     peak1 = min(find(trajdata1{maxtraj}(:,5) == out.peakrate1));
     peak2 = min(find(trajdata2{maxtraj}(:,5) == out.peakrate2));
     
     out.peakdist = (abs(peak2-peak1))*binsize; % in cm

 else  % Different max traj
     
     include = [];
     for i=1:length(maxtraj)
         currtraj = maxtraj(i);
         
         peakrate1(i) = nanmax(trajdata1{currtraj}(:,5));
         peakrate2(i) = nanmax(trajdata2{currtraj}(:,5));
         
         peak1 = min(find(trajdata1{currtraj}(:,5) == peakrate1(i)));
         peak2 = min(find(trajdata2{currtraj}(:,5) == peakrate2(i)));
         peakdist(i) = (abs(peak2-peak1))*binsize;
         
         if (peakrate1(i)>=thresh && peakrate2(i)>=thresh && ~isnan(peakdist(i))...
                 && length(find(isnan(trajdata2{currtraj}(:,5))))<=0.5*length(trajdata2{currtraj}(:,5))...
                 && length(find(isnan(trajdata1{currtraj}(:,5))))<=0.5*length(trajdata1{currtraj}(:,5)));
             include = [include,i];
         end
         
     end
     
     
     if length(include)==1 % if only 1 traj is above thresh
        out.peakrate1 = peakrate1(include);
        out.peakrate2 = peakrate2(include);
        out.peakdist = peakdist(include);
     end       
         
     if isempty(include) || length(include)>1 % if neither or both is greater than thresh, take the mean
        out.peakdist = nanmean(peakdist);
        out.peakrate1 = nanmax(peakrate1);
        out.peakrate2 = nanmax(peakrate2);
     end
     
     
 end
 






