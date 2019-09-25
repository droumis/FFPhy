function out = calcaveragetrace(index, excludetimes, eeg,ripple, varargin)
% function out = calcaveragetrace(index, excludetimes, ripple, varargin)
%
% set options
appendindex = 1;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% assign a temporary variable for eeg
e = eeg{index(1)}{index(2)}{index(3)};
r = ripple{index(1)}{index(2)}{index(3)};
times = geteegtimes(e);

clear eeg;
clear ripple;

% apply excludetimes
includetimes = getincludedtimes(excludetimes);

%Find the start and end indices for valid times
startind = lookup(includetimes(:,1),times)-150;
endind   = lookup(includetimes(:,2),times)+150;

if ~isempty(startind)
    while startind(1)<0
        startind = startind(2:end);
        endind = endind(2:end); 
    end
    while endind(end)>size(r.data,1)
        startind = startind(1:end-1);
        endind = endind(1:end-1);
    end
end

s = zeros(length(startind),2);
% find index of peak of ripple power
for i = 1:length(startind)
    tmprip = double(r.data(startind(i):endind(i),3));
    [maxtab mintab] = peakdet(tmprip,0.25);
    peak = maxtab(find(maxtab(:,2)==max(maxtab(:,2))),1);
    s(i,1) = startind(i)+peak(1) - 250;
    s(i,2) = startind(i)+peak(1) + 250;
end

% get average trace during goodtimes
rip = zeros(size(s,1),501);
eeg = zeros(size(s,1),501);
for i = 1:size(s,1)
    tmpeeg = e.data(s(i,1):s(i,2),1);
    eeg(i,:) = tmpeeg;
    tmprip = double(r.data(s(i,1):s(i,2),3));
    rip(i,:) = tmprip;
end

if size(eeg,1)>1
    eeg = mean(eeg);
    rip = mean(rip);
end

out.eeg = eeg;
out.rip = rip;
out.est = (eeg(rip==max(rip))-mean(e.data))./std(e.data);

if (appendindex)
    out.index = index;
end
