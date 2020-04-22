function f = runfilter(f, varargin)

useparpool = 0; %DR added flag to activate parallel pool for multi animal speedup
outputDayEpTetCells = 0;
if ~isempty(varargin)
    assign(varargin{:});
end

if useparpool
    parfor an = 1:length(f) %this will negate any nested parfor loops
        iterator = f(an).iterator;
        f(an) = feval(iterator,f(an), 'outputDayEpTetCells', outputDayEpTetCells);
    end
else    
    for an = 1:length(f)
        iterator = f(an).iterator;
        f(an) = feval(iterator,f(an), 'outputDayEpTetCells', outputDayEpTetCells);
    end
end