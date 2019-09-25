function indices= findIndex(t, otimes, lo_or_hi)

% find timestep index corresponding to time t in _ordered_ list otimes
% uses bisection search
% lo_or_hi= 1 (lo), 2 (hi)


if nargin < 3
    lo_or_hi= 1; % = lo
else
    if lo_or_hi~= 1 & lo_or_hi ~= 2
	error('incorrect choice');
    end
end

nt= length(t);
no= length(otimes);
indices= zeros(nt,1);
b=[1, no];
for i=1:nt
    [b(1) b(2)]= auxrun(t(i), otimes, 1, no);
    indices(i)= b(lo_or_hi);
end


function [lo, hi]= auxrun(t, otimes, lo_ini, hi_ini)

small= (otimes(2)-otimes(1))/1000;
lo= lo_ini; hi= hi_ini;

if (t < otimes(lo)-small) | (t > otimes(hi)+small)
    str= sprintf('%f, (%f, %f)', t, otimes(lo), otimes(hi));
    error(['requested time out of bounds: ' str]);
end
if t== otimes(lo)
    hi= lo;
    return;
elseif t== otimes(hi)
    lo= hi;
    return
end


while hi-lo > 1
    mid= floor((hi+lo)/2);
    if t == otimes(mid)
	lo= mid;
	hi= mid;
	break;
    elseif t > otimes(mid)
	lo= mid;
    else 
	hi= mid;
    end
end
