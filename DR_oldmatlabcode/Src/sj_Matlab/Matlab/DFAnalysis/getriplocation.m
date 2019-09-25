function [out] = getriplocation(index, excludeperiods, ripples, pos, linpos, includestate, varargin)
% out = getriplocation(index, excludeperiods, ripples, pos, linpos,  options)
%
% reports the position for all times when ripples detected
%
%   options are
%	'minenergy', E
%		     specifies the minimum energy of a valid ripple event
%   'numtetrodes'
%           specifies number of tetrodes a ripple must be recorded on to be
%           included in analysis, default 1
%   'proptetrodes', examples: 1, 0.5, 0.25,
%           proportion of tetrodes a ripple must be recorded on to be
%           included in analysis
%   'appendindex' , 1 or 0, default 0
%           set to 1 to append the day and epoch index to the output [day
%           epoch]
%
%   out is out.incllindist traj & lindist included by included times: [ traj lindist]
%          out.riplindist traj & lindist of ripple times: [traj lindist]


% assign the options
numtetrodes = 1;
minenergy = 0;
proptetrodes = [];
appendindex = 0;
plotfig = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'numtetrodes'
            numtetrodes = varargin{option+1};
        case 'proptetrodes'
            proptetrodes = varargin{option+1};
        case 'appendindex'
            appendindex = varargin{option+1};
        case 'minenergy'
            minenergy = varargin{option+1};
        case 'plotfig'
            plotfig = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if ~isempty(proptetrodes)
    numtetrodes = round(size(index,1) * proptetrodes); %size(index,1) is total number of tetrodes input
end

tetlist = index(:,3);
r = ripples{index(1,1)}{index(1,2)}{tetlist(1)};

postimes = pos{index(1,1)}{index(1,2)}.data(:, 1);
times = postimes;
%times = r.timerange(1):0.001:r.timerange(end);% sampled every millisecond
nrip = zeros(size(times));
for t = 1:length(tetlist)
    tmprip = ripples{index(1,1)}{index(1,2)}{tetlist(t)};
    if ~isempty(tmprip.starttime)  %is any ripples
        if (minenergy == 0)
            % get all the times
            rtimes = [tmprip.starttime tmprip.endtime];
        else
            % get the indeces for the ripples with energy above minenergy
            rvalid = find(tmprip.energy > minenergy);
            rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
        end
        % create another parallel vector with bordering times for zeros
        nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
        rtimes = reshape(rtimes', length(rtimes(:)), 1);
        rtimes(:,2) = 1;
        nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
            length(nrtimes(:)), 1) ; r.timerange(2)];
        nrtimes(:,2) = 0;
        % create a new list with all of the times in it
        tlist = sortrows([rtimes ; nrtimes]);
        % use interp to create a set of ones and zeros for each time
        % and add to nrip to get a cumulative count of the number of
        % ripples per timestep
        try
            nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest'); %number of ripples on each tetrode
        catch
            keyboard
        end
    end
end

%apply excludetimes to nrip
includetimes = ~isExcluded(times', excludeperiods); %list of ones and zeros
if size(nrip) ~= size(includetimes)
    includetimes = includetimes';
end
nrip(includetimes==0) = 0; %exclude ripples if excluded by excludetimes
inclripindex = (nrip >= numtetrodes); %list of 0 and 1s, 1 indicate times with ripple that meets criteria, sampled as postimes

% %find times with ripples and change to pos sampling
% inclriptimesmsec = (nrip >= numtetrodes); %list of 0 and 1s, 1 indicate times with ripple that meets criteria, sampled every msec
% ripexcludeperiods = getExcludePeriods(times', inclriptimesmsec');
% postimes = pos{index(1,1)}{index(1,2)}.data(:, 1);
% inclripindex = ~isExcluded(postimes, ripexcludeperiods); %list of 0 and 1s, 1 indicate times with ripple that meets criteria, sampled in pos time

%get rip location at each time that is included
allposx = pos{index(1,1)}{index(1,2)}.data(:, 2); % x position for all times
allposy = pos{index(1,1)}{index(1,2)}.data(:, 3); % y position for all times
%includetimes = ~isExcluded(postimes, excludeperiods);
inclposx = pos{index(1,1)}{index(1,2)}.data(includetimes, 2); % x position for times not excluded by excludeperiods
inclposy = pos{index(1,1)}{index(1,2)}.data(includetimes, 3); % y position for times not excluded by excludeperiods
ripposx = allposx(inclripindex, 1); % x position for ripple times
ripposy = allposy(inclripindex, 1); % y position for ripple times
[state, lindist] = getbehavestate(linpos, index(1,1), index(1,2), includestate, 'minlinvelocity', 0);
alltraj = state;
alllindist = lindist;
incltraj = alltraj(includetimes, :); % traj for times not excluded by excludeperiods
incllindist = alllindist(includetimes, :); % lindist for times not excluded by excludeperiods
riptraj = alltraj(inclripindex, 1);
riplindist = alllindist(inclripindex, 1);

%for debugging.  fraction should be less than 1
maxlind = 165; %max(alllindist(ismember(alltraj, 1:6)));
edges = 0:20:maxlind+20;
for t = 1:6
    inclloclin = incllindist(incltraj == t);
    riploclin = riplindist(riptraj==t);
    if any(inclloclin) && any(riploclin)
        h = histc(inclloclin,edges);
        riph = histc(riploclin,edges);
        if size(h,1)==1
            h = h';
        end
        if size(riph,1)==1
            riph = riph';
        end
        if rem(t,2) == 0
            h(1:end-1) = flipud(h(1:end-1));
            riph(1:end-1) = flipud(riph(1:end-1));
        end
        if any((riph./h)>1)
            error('fraction>1')
        end
    end
end

% plot
if plotfig == 1
    figure(1)
    hold on
    plot(allposx, allposy, 'color', [0.5  0.5 0.5])
    plot(inclposx, inclposy, '.b', 'markersize', 15)
    plot(ripposx, ripposy, '.r', 'markersize', 15)
    title(['D', num2str(index(1,1)), ' E', num2str(index(1,2)), 'grey: allpos, blue:includedpos, red:rippos'])
    
    maxlind = 165; %max(alllindist(ismember(alltraj, 1:6)));
    edges = 0:20:maxlind+20;
    for t = 1:6
        inclloclin = incllindist(incltraj == t);
        riploclin = riplindist(riptraj==t);
        if any(inclloclin)
            h = histc(inclloclin,edges);
            riph = histc(riploclin,edges);
            if size(h,1)==1
                h = h';
            end
            if size(riph,1)==1
                riph = riph';
            end
            if rem(t,2) == 0
                h(1:end-1) = flipud(h(1:end-1));
                riph(1:end-1) = flipud(riph(1:end-1));
            end
            figure(2)
            subplot(6,1,t)
            bar(edges+mean(diff(edges))/2, h./sum(h), 'k')
            hold on
            if ~isempty(riph)
                bar(edges+mean(diff(edges))/4, riph./sum(riph),'r')
                hold on
            end
            title(['Traj', num2str(t)])
            
            if ~isempty(riph)
                figure(3)
                subplot(6,1,t)
                plot(edges+mean(diff(edges))/4, riph./h,'b', 'linewidth',3)
                hold on
                title(['Traj', num2str(t)])
                
                if any((riph./h)>1)
                    error('fraction>1')
                end
            end
            
        end
    end
    
    figure(2)
    xlabel('Distance (cm)')
    ylabel('Fraction')
    subtitle(['D', num2str(index(1,1)), ' E', num2str(index(1,2)), ' black:includedpos, red:rippos'])
    
    figure(3)
    xlabel('Distance (cm)')
    ylabel('Fraction of time spent in SWR')
    subtitle(['D', num2str(index(1,1)), ' E', num2str(index(1,2)), ])
    
    pause
    clf(1)
    clf(2)
    clf(3)
end

%output
if appendindex ==1
    appendcol = [index(1,1:2)];
else
    appendcol = [];
end
out.incllindist = [repmat(appendcol, length(incltraj), 1) incltraj incllindist]; %traj & lindist included by included times: [ day epoch traj lindist]
out.riplindist =  [repmat(appendcol, length(riptraj), 1) riptraj riplindist]; %traj & lindist of ripple times: [ day epoch traj lindist]




end

