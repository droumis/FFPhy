function outnum = roundn(innum,varargin)
% Round numbers to a specified power of 10.
% Syntax
% outnum = roundn(innum)
% outnum = roundn(innum,n)
% Description
% outnum = roundn(innum) rounds the elements of innum to the nearest one-hundredth.
% outnum = roundn(innum,n) specifies the power of 10 to which the elements of innum are rounded. For example, if n = 2, round to the nearest hundred (10^2).

numvarargs = length(varargin);
if numvarargs > 1
    error('roundn','requires at most 1 optional input');
end

optargs = {-2};

% skip any new inputs if they are empty
newVals = ~cellfun('isempty', varargin);
optargs(newVals) = varargin(newVals);

% Place optional args in named variables
[n] = optargs{:};
factor = 10^(-n);
outnum = round(innum*factor)/factor;