function out = stack(inputA, inputB)
%out = stack(inputA, inputB)
%Creates a matrix with the two matrices stacked.  If the 2nd
%dimension of the two matrices are equal, this is equivalent to
%[inputA;inputB].  If not, then the smaller matrix is filled in with NaN's
%until the sizes match.

if ~(isnumeric(inputA) & isnumeric(inputB))
    error('Both inputs must be numeric');
end

sizA = size(inputA,2);
sizB = size(inputB,2);
try
    out = [inputA;inputB];
catch

    if (sizA > sizB)
        inputB(:,sizB+1:sizA) = nan;
        out = [inputA;inputB];
    elseif (sizB > sizA)
        inputA(:,sizA+1:sizB) = nan;
        out = [inputA;inputB];
    else
        error('Can not stack');
    end
end

