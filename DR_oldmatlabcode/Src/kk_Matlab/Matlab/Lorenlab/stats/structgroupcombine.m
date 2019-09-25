function out = structgroupcombine(f)
%out = structgroupcombine(f)
%Combines data for each group across animals.  It is
%assumed that the struct data is stored in cells in the filter's output field.

   
out = [];
for an = 1:length(f)
    for g = 1:length(f(an).output)
        
        tmp = f(an).output{g};
        if ((length(out) < g)|(isempty(out)))
            out{g} = tmp;
        else
            out{g} = [out{g} f(an).output{g}];
        end
        
    end
end