function [out] = sj_getandsave_ripacrosstet(prefix,day)
% Shantanu - From DFAsj_getriprate. This is done outside the filter framework
% Do this once - Load "ripples" file for day, eg. REdripples03, and extract
% ripple detection/startime/rate/etc across all riptetlist.
% Save it in a "ripplesall" file, eg REdripplesall03 - Aug2012. Make it "ripplall"
% out = sj_getandsave_ripacrosstet('REd',1);

% Use only 3 SD
% Criterion of 1 tet. save as ripplesalltet{day}{ep}{1}
% Criterion of 2 tet. save as ripplesalltet{day}{ep}{2}

%   'minthresh', minthresh
%		     specifies a minimum threshold in stdev units for a valid
%			ripple event  (default 0)
%   'numtetrodes'
%           specifies number of tetrodes a ripple must be recorded on to be
%           included in analysis, default 1

% NEED TO ADD A STIMFILTER

% -------------------------------------------
% Set directories
% -------------------------------------------

switch prefix
    case 'HPa'
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
        dire = '/data25/sjadhav/HPExpt/HPa';
        riptetlist=[1,4,5,6];
      
    case 'REc'
        directoryname = '/data25/sjadhav/RippleInterruption/REc_direct/';
        dir = '/data25/sjadhav/RippleInterruption/REc';
        riptetlist = [4,6,8,9,10];
    case 'REd'
        directoryname = '/data25/sjadhav/RippleInterruption/REd_direct/';
        dir = '/data25/sjadhav/RippleInterruption/REd';
        riptetlist = [3,4,5,6,10,11];
    case 'REe'
        directoryname = '/data25/sjadhav/RippleInterruption/REe_direct/';
        dir = '/data25/sjadhav/RippleInterruption/REe';
        riptetlist = [3,4,6,11,12,13];
    case 'REf'
        directoryname = '/data25/sjadhav/RippleInterruption/REf_direct/';
        dir = '/data25/sjadhav/RippleInterruption/REf';
        riptetlist = [1,5,9,10,11,12];
    case 'RCa'
        directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct/';
        dir = '/data25/sjadhav/RippleInterruption/RCa';
        riptetlist = [2,3,4,6,9];
    case 'RCb'
        directoryname = '/data25/sjadhav/RippleInterruption/RCb_direct/';
        dir = '/data25/sjadhav/RippleInterruption/RCb';
         riptetlist = [3,4,9,10,11,12];
    case 'RCc'
        directoryname = '/data25/sjadhav/RippleInterruption/RCc_direct/';
        dir = '/data25/sjadhav/RippleInterruption/RCc';
        riptetlist = [3,4,5,6,11,13];
    case 'RCd'
        directoryname = '/data25/sjadhav/RippleInterruption/RCd_direct/';
        dir = '/data25/sjadhav/RippleInterruption/RCd';
        riptetlist = [1,2,3,4,5,6];
end

tetlist = riptetlist;
% Get rippleep1 or ripple file
ripfile = sprintf('%s/%sripplesep1%02d.mat', directoryname, prefix, day);
load(ripfile);

% Filter for stim times if necessary - for Rip Disruption Expt
% Tet filter for ripple detection
riptetfilter = '(isequal($descrip, ''riptet''))';
%timefilter = {{'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter}};
f(1).animal{1} = prefix;
f(1).animal{3} = prefix;
f(1).animal{2} = directoryname;

for ep = 1:length(ripples{day})
    ep;
    f(1).epochs{1} = [day ep];
    %f = setfiltertime(f, timefilter);
    %excludeperiods = f(1).excludetime{1}{1};
    excludeperiods=[];
    
    r=ripples{day}{ep}{riptetlist(1)}; % Get time base from 1st tet
    
    times = r.timerange(1):0.001:r.timerange(end);
    nrip = zeros(size(times));
    nrip_size = zeros(length(tetlist),length(times));
    
    for t = 1:length(tetlist)
        tmprip = ripples{day}{ep}{tetlist(t)};
        
        % get all the times
        rtimes = [tmprip.starttime tmprip.endtime];
        rvalid = 1:length(tmprip.maxthresh); % All ripples are valid
        
        % For ripsize
        currsize=tmprip.maxthresh(rvalid); % Size of each ripple on current tetrode
        
        % create another parallel vector with bordering times for zeros
        nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
        rtimes = reshape(rtimes', length(rtimes(:)), 1);
        rtimes(:,2) = 1;
        
        % For ripsize
        % Make currsize match the format of rtimes and add it there
        currsize(:,2) = currsize;
        currsize = reshape(currsize', length(currsize(:)), 1);
        rtimes(:,3) = currsize;
        %     rtimes_size = rtimes(:,1);
        %     rtimes_size(:,2) = currsize;
        
        % Update nrtimes
        nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
            length(nrtimes(:)), 1) ; r.timerange(2)];
        nrtimes(:,2) = 0;
        nrtimes(:,3) = 0;
        
        % Create a new list with all of the times in it with 1 or 0 for ripple
        % or not
        tlist = sortrows([rtimes ; nrtimes]);
        % Alternate method
        %     % create a new list with all of the times in it with ripsize instead of
        %     % 1 or 0
        %     tlist_size = sortrows([rtimes_size ; nrtimes]);
        
        % use interp to create a set of ones and zeros for each time
        % and add to nrip to get a cumulative count of the number of
        % ripples per timestep
        try
            nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
        catch
            keyboard
        end
        
        try
            nrip_size(t,:) = interp1(tlist(:,1), tlist(:,3), times, 'nearest');
        catch
            keyboard
        end
        % mean, std and threshold of ripple envelope used for each tet in current epoch
        % Return a mean acroos all tets for comparing acroos tets later
        % Also, sj_calcriprate can figure out which tet is assigning ripsize,
        % and return its baseline, std and threshold. % Not implented now
        rip_baseline(t) = tmprip.baseline;
        rip_std(t) = tmprip.std;
        rip_thresh(t) = tmprip.threshold;
    end % end tetlist
    
    % calculate ripple rate, etc\
    rip.times = times; %sampled every millisecond
    clear times;
    rip.nripples = nrip; %number of ripples on all tetrode
    clear nrip;
    rip.ripsize = nrip_size;
    clear nrip_size;
    for numtetrodes = 1:2
        if numtetrodes==1
            ripplesall{day}{ep}{numtetrodes}.baseline = mean(rip_baseline);
            ripplesall{day}{ep}{numtetrodes}.std = mean(rip_std);
            ripplesall{day}{ep}{numtetrodes}.threshold = mean(rip_thresh);
            ripplesall{day}{ep}{numtetrodes}.times = rip.times; %sampled every millisecond
            ripplesall{day}{ep}{numtetrodes}.nripplesvec = rip.nripples; %number of ripples per ms combined on all tetrodes
        end
        out = sj_calcriprate(rip, excludeperiods, numtetrodes);
        ripplesall{day}{ep}{numtetrodes}.starttime = out.ripsttime;
        ripplesall{day}{ep}{numtetrodes}.maxthresh = out.ripsize;
        ripplesall{day}{ep}{numtetrodes}.ripntet = out.ripntet;
        ripplesall{day}{ep}{numtetrodes}.riprate = out.rip(1);
        ripplesall{day}{ep}{numtetrodes}.rippercenttime = out.rip(2);
        ripplesall{day}{ep}{numtetrodes}.stilltime = out.rip(3);
        ripplesall{day}{ep}{numtetrodes}.indivthresh = out.ripsize_indiv;
    end % end numtetrodes
    clear rip_baseline rip_std rip_thresh
 
end % end ep
save(sprintf('%s%sripplall%02d.mat', directoryname, prefix, day), 'ripplesall');



