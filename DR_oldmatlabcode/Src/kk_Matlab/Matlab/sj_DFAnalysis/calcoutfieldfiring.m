function out = calcoutfieldfiring(index, excludeperiods, spikes, linpos, varargin)
%Calculates out-of-field firing, based on the linearized occ-normalized
%firing rates.  Out-of-field is defined as bins with rates between .05 and
%1 Hz.

appendindex = 0;
thresh = 0;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'thresh'
                thresh = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end
meanrate = calctotalmeanrate(index, excludeperiods, spikes);
if (meanrate >= thresh)
    linfields = filtercalclinfields(index, excludeperiods, spikes, linpos);
    rates = [];
    for i = 1:length(linfields.trajdata)
        rates = [rates; linfields.trajdata{i}(find(~isnan(linfields.trajdata{i}(:,5))),5)];
    end
    numoofbins = length(find((rates > .05) & (rates < 1) ));
    out = numoofbins/length(rates);
else
    out = nan;
end