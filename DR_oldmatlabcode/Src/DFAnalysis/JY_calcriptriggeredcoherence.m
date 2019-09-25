
function out = JY_calcriptriggeredcoherence(index, excludetimes, data, eeg, ripples,varargin)
% function out = calccoherence_test(index, excludetimes, eeg, varargin)
%  Based on Maggie's calccoherence.m
%  Plots the coherence for an eeg tetrode pair. If you use a time filter,
%  excluded times are removed and the includedtimes are averaged together.
%
%   out is a structure with the following fields
%       coherence-- This is the coherence for the tetrode pairs
%       frequency-- Frequency vector
%       index-- Only if appendindex is set to 1 (default)



for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'Fs'
                params.Fs = varargin{option+1};
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'err'
                params.err = varargin{option+1};
            case 'rippletetrode'
                ripref=varargin{option+1};
                
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

% get ripple information
minstd=5;
getripindex=find(ripples{index(1)}{index(2)}{ripref}.maxthresh>=minstd);
%rippletimes=ripples{index(1)}{index(2)}{ripref}.midtime(getripindex)*10000;
ripplesize=ripples{index(1)}{index(2)}{ripref}.maxthresh(getripindex);
rippletimes=ripples{index(1)}{index(2)}{ripref}.starttime(getripindex)*10000;


rippletimes=rippletimes(~isExcluded(rippletimes, excludetimes));


% assign a temporary variable for eeg
e1 = eeg{index(1)}{index(2)}{index(3)};
e2 = eeg{index(1)}{index(2)}{index(4)};
e1start = e1.starttime;
e2start = e2.starttime;

e1times = geteegtimes(e1);
e2times = geteegtimes(e2);

e1sampfreq=e1.samprate;

if length(e1times)>length(e2times)
    temp = lookup(e2times,e1times);
    e1 = e1.data(temp);
    e2 = e2.data;
elseif length(e2times)>length(e1times)
    temp = lookup(e1times,e2times);
    e1 = e1.data;
    e2 = e2.data(temp);
elseif length(e1times)==length(e2times)
    e1 = e1.data;
    e2 = e2.data;
end

% put window around ripples for analysis


halfwindow=10000;
i=1;
e1trials=[];
e2trials=[];
while i<length(rippletimes)
    % go through each tetrode
    tstart=rippletimes(i,1)-halfwindow;
    tend=rippletimes(i,1)+halfwindow;
    %tend=rippletimes(i,2);+halfwindow;
    
    %ripple data
    
    
    tstartindex=lookup(tstart,e1times*10000);
    if tstartindex ==1 % check for mismatch between epoch start times defined by times.mat and position reconstruction
       i=i+1;
    end
    tendindex=lookup(tend,e1times*10000);
    
    if tendindex==size(e1,1);
       i=size(rippletimes,1)+1;
    end
    
    e1trials(i,:)=e1(tstartindex:tendindex);
    e2trials(i,:)=e2(tstartindex:tendindex);
    i=i+1;
end

% parse the options
params = [];
params.Fs = 1500;
params.fpass = [1 300];
params.err = [1 0.05];
appendindex = 1;
window =[0.5 0.005];
params.tapers=[3 5];


% compute full coherence
[C,phi,S12,S1,S2,times,freq]=cohgramc(e1trials',e2trials',window, params);
% 
% figure;
% meanC=mean(C,3);
% imagesc(t-halfwindow/10000,f,meanC');
% axis xy;
% colorbar;

%apply excludetimes
% goodtimes = ~isExcluded(t+e1start, excludetimes);
% tempcohere = C;
% tempcohere= tempcohere(goodtimes,:);
% coherence = tempcohere;%mean(tempcohere);

p = data{index(1)}{index(2)}.Pos.correcteddata;
time = lookup(times+e1start, p(:,1));
speed = p(time,5);

out.times = times;
out.coherence = C;
out.frequency = freq;

if (appendindex)
    out.index = index;
end

end
