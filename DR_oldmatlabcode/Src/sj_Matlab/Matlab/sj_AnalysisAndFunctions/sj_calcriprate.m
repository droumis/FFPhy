function out = sj_calcriprate(rip, excludeperiods, numtetrodes, varargin)
%out = calcriprate(rip, excludetimes, numtetrodes)
%
%   rip.times: time of each sample in nripples, 1msec bins
%   rip.nripples: number of electrodes with a ripple recorded at that time
%   point
%   rip.ripsize: Ripple Size per ms - same size is reported from start to end of ripple   
%
%   Not Implemented Yet:
%   sj_calcriprate can figure out which tet is assigning ripsize, 
%   and return its baseline, std and threshold.  
%   Baseline, Std and Threshold for eac tet in epoch used to calculate nrip
%   and nsize: rip.baseline, rip.std, rip.thrs 
%   
%
%   out is [ rate proportiontime] 
%       proprotiontime is proportion of included time during which ripples
%   were recorded
%       rate is number ripples/sec during included time
%
% If DIO stimtimes are given, exclude them from given time ranges

stimtime = [];
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'stimtime'
            stimtime = varargin{option+1};
    end
end

%apply excludetimes to nripples 
includetimes = ~isExcluded(rip.times, excludeperiods); %list of ones and zeros sampled every millisecond, ones = included, zeros = excluded
if size(rip.nripples) ~= size(includetimes)
    includetimes = includetimes';
end
includerips = rip.nripples .* includetimes;
if isfield(rip,'ripsize')   
    %includeripsize = max(rip.ripsize,[],1) .* includetimes; % Take maximum ripple size across tetrodes
    includeripsize_indiv = rip.ripsize .* repmat(includetimes,size(rip.ripsize,1),1);
end

%calculate proportion of time spent in ripples satisfying ntet condition
riptimecount = length( find(includerips >=numtetrodes) );
%rippercenttime = riptimecount / sum(includetimes); %proportion total recorded time that included ripples 
% Aug 2012 - Change. THis should be fraction of all times, not just includetimes.
rippercenttime = riptimecount / length(rip.times); %proportion total recorded time that included ripples 
% Also get [total time - excludetimes]
stilltime = sum(includetimes); % This has usually already been filtered for speed for ripples 

%calculate number of ripples per time 
a = zeros(1, length(includerips));
a(find(includerips >= numtetrodes)) = 1;
ripcount = length(find(diff(a) == 1));

% Calculations for ntet and ripsize 
% Get ntet for each ripple
ripst = find(diff(a)==1) + 1;
ripend = find(diff(a)==-1) - 1;
if length(ripst)~=length(ripend)       % If lengths dont match
    if length(ripst) > length(ripend)  % Ripple at end of vector a
        ripend = [ripend, length(a)];
    else                               % Ripple at start of vector a
        ripst = [1, ripst];
    end
end

% If lengths still dont match, 1) you can trim vectors to match lth
% 2) Alternatively, use try-catch to catch error
% 1)
% if length(ripst)~=length(ripend)
%     comlth = min([length(ripst) length(ripend)]);
%     ripst = ripst(comlth);
%     ripend = ripend (comlth)
% end
%ripstend = [ripst; ripend];

%2)
ripsttime = rip.times(ripst); % Ripple start times
ripntet = zeros(size(ripst));
ripsize = zeros(size(ripst));
ripsize_indiv = zeros(size(ripst,2),size(rip.ripsize,1));
try
    for ri=1:length(ripst),
        if ripend(ri) > ripst(ri)
            ripntet(ri) = max(includerips(ripst(ri):ripend(ri))); % Max no of tetrodes on which ripple was detected 
        else
            ripntet(ri) = includerips(ripst(ri));
        end
        
        % For size: just use start of ripple - same size reported
        % throughout ripple
        if isfield(rip,'ripsize')  
            %ripsize(ri) = includeripsize(ripst(ri));
            ripsize_indiv(ri,:) = includeripsize_indiv(:,ripst(ri));
            ripsize(ri) = max(includeripsize_indiv(:,ripst(ri)));
        end
    end
catch
    keyboard
end
%ripntet = includerips(ripidx+1);
%length(find(ripntet>=2))

% if DIO/stim exists, apply exclude periods to it and remove it
if isempty(stimtime)
    riprate = ripcount / (sum(includetimes)/1000);  %divide sum of all included time by 1000 change from msec to seconds
    % Aug 2012 - Change. THis should be fraction of all times, not just includetimes.
    riprate = ripcount / (length(rip.times)/1000);  %divide sum of all time by 1000 change from msec to seconds
else
    
    includetimepoints = ~sj_isTimePtExcluded(stimtime./1000, excludeperiods);
    if size(stimtime) ~= size(includetimepoints)
        includetimepoints = includetimepoints';
    end
    stimtime = stimtime(includetimepoints);
    
    % Find stimulations that cause a artifactual ripple
    sind = lookup(stimtime/1000,ripsttime); % Find the closest ripple start times
    tdiff = stimtime' - ripsttime(sind)*1000; % See how far ripple is from stimulation in ms
    removestimidx = sind(find(abs(tdiff)<=100)); % If less than 100ms, artifactual ripple has to go
    
    % First adjust ripcount
    ripcount = ripcount - length(removestimidx);
    riprate = ripcount / (sum(includetimes)/1000);  %divide sum of all included time by 1000 change from msec to seconds
    % Aug 2012 - Change. THis should be fraction of all times, not just includetimes.
    riprate = ripcount / (length(rip.times)/1000);  %divide sum of all time by 1000 change from msec to seconds

    % Then adjust ripntet and ripsize 
    ripsize(removestimidx) = [];
    ripntet(removestimidx) = [];
    ripsize_indiv(removestimidx,:) = [];
end
    
    
%out = [riprate rippercenttime*100 stilltime/1000]; 
out.rip = [riprate rippercenttime*100 stilltime/1000];   % rippercentime in percent, stilltime in seconds
out.ripntet = ripntet;
out.ripsize = ripsize; % Vector of zeros if size field does not exist
out.ripsttime = ripsttime;
out.ripsize_indiv = ripsize_indiv;
