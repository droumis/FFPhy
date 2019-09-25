function [ind, tfirst]= findBursts(t, timesep)
% function [indices]= getBursts(spiketimes)
% from a given list of spike times select the first spikes in each burst
%   timesep= 0.2    [sec] minimum separation between bursts
% how?

if nargin < 1 | isempty(t); indices= zeros(0,1); return; end
if nargin < 2; timesep= 0.2; end

if size(t,1) == 1 
    id=2;
elseif  size(t,2) == 1
    id=1;
else
    error('cannot deal with matrices');
end


d= diff(t);
if id== 1
    ind= [0;find(d > timesep)];
else
    ind= [0,find(d > timesep)];
end

ind= ind+1;
if nargout > 1; tfirst= t(ind); end
