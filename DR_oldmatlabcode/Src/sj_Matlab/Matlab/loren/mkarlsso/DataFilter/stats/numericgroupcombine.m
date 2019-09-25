function out = numericgroupcombine(f, appendindex)
%out = numericgroupcombine(f)
%out = numericgroupcombine(f,appendindex)
%Combines data for each group across animals.  It is
%assumed that the numeric data is stored in cells in the filter's output field.
%If the 2nd dimension does not match, then the difference is filled in with
%NaN's. If appendindex is 1, the animal number is appended to each entry
%(default 1).

if (nargin < 2)
    appendindex = 1;
end
    
out = [];
for an = 1:length(f)
    for g = 1:length(f(an).output)
        if (length(out) < g)
            out{g} = [];
        end
        tmp = f(an).output{g};
        if (appendindex)
            firstcolumn = ones(size(f(an).output{g},1),1) * an;
            out{g} = stack(out{g}, [firstcolumn f(an).output{g}]);
        else
            out{g} = stack(out{g}, f(an).output{g});
        end
    end
end
   