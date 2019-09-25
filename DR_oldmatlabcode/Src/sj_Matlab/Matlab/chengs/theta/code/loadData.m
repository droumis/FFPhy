function [data,T]= loadData(num, bursts)
% function [data,T]= loadData(num)

if nargin < 2 | isempty(bursts); bursts= 0; end
global fmaux behavdata spikedata

d=num(1); e= num(2); t= num(3); c= num(4);

loadVar(fmaux.data2dir, 'behavdata', d);

if ~isfield(fmaux, 'datafile') | isempty(fmaux.datafile)
    loadVar(fmaux.data2dir, 'spikedata', d);
else
    loadVar(fmaux.data2dir, fmaux.datafile, d);
end


data= behavdata{d}{e};
data.spiketimes= spikedata{d}{e}{t}{c}.time;
data.spikeindex= spikedata{d}{e}{t}{c}.index;
if bursts
    ibursts= findBursts(data.spiketimes);
    data.spiketimes= data.spiketimes(ibursts);
    data.spikeindex= data.spikeindex(ibursts);
end

if nargout > 1
    T= length(data.time);
end
