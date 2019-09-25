function out = calcgroupstats(cellarray,columns)
% out = calcgroupstats(cellarray,columns)
% Computes simple stats of each column of each matrix inside cellarray{1:N}
% Output is a structure of length N

out = [];
for i = 1:length(cellarray) 
    if isnumeric(cellarray{i})      
            if (nargin < 2)
                columns = 1:size(cellarray{i},2);
            end
            out(i).mean = noNanMean(cellarray{i}(:,columns));         
            out(i).stderr = noNanStderr(cellarray{i}(:,columns));         
            out(i).bootstderr = noNanBootStderr(cellarray{i}(:,columns));
            out(i).groupnum = i;
    else
        error('The contents of each cell must be numeric');
    end
end