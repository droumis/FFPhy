function po = bootstrap_smoothphase(num_bootstraps,p,x,xi,varargin);
%
%   po = bootstrap_smoothphase(num_bootstraps,p,x,xi,[params])
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERROR-CHECKING and
% FUNCTION-SPECIFIC CONSTANTS 

if length(varargin)==1
    params = varargin{1};
elseif length(varargin)>1
    error('too many extra arguments!')
end

for i = 1:num_bootstraps
    idx = unidrnd(numel(x),size(x));
    if exist('params')
        po(:,i) = smoothphase(p(idx),x(idx),xi,params);
    else
        po(:,i) = smoothphase(p(idx),x(idx),xi);
    end
end

