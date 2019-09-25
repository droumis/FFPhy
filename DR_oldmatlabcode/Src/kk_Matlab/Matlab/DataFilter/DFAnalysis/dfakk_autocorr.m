function out = dfakk_autocorr(ind,eind,excludetimes, spikes, linpos, thetagnd, varargin)

% kk 7.16.13, function is derived from dfakk_calcxcorrmeasures

%function out = calcxcorrmeasures(index, excludetimes, spikes, varargin)
% Calculates the excess correlation and RMS time lag for the specified cell
% pairs using only spikes not excluded by the excludetimes
% 
%
% Options:
%   'bin', n 	binsize in sec. Default 0.001 (1 ms)
%   'tmax', n	maximum time for cross correlation. Default 1 sec.


bin = 0.001;
tmax = 1;       


for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'bin'
                bin = varargin{option+1};
            case 'tmax'
                tmax = varargin{option+1};    
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end


% set the defaults
out.time = [];
out.ac1 = [];
out.index = ind;



% for each cell we calculate the cross correlation 
try
    t1 = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data;
catch
    % if either of those produced an error, we return NaNs 
    return;
end

% exclude spikes in excluded periods
t1inc = [];
if (length(t1))
    t1inc = t1(find(~isExcluded(t1(:,1), excludetimes)),1);
else
    return
end



% autocorrelation
ac1 = spikexcorr(t1inc, t1inc, bin, tmax);
out.ac1 = ac1.c1vsc2;
out.time = ac1.time;


