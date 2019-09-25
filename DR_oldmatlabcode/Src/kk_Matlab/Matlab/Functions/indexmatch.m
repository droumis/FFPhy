function out = indexmatch(M1,M2,indexcolumns,varargin)
% out = indexmatch(M1,M2,indexcolumns)
% M1 and M2 are matrices where the columns INDEXCOLUMNS represent an
% index identifier to the data in each row.  For any rows with matching
% indices, the functions add a row to the output with the index and the
% remaining data from M1 and M2.  Warning: it is assumed that there will be
% either 0 or 1 matches for each index.  If there are more than one, it
% will pick the first one from M2.


outindex = [];
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'outindex'
                outindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end


indexcolumns = (indexcolumns(:))'; %make sure indexcolumns is a horizontal vector
M1datacol = [max(indexcolumns)+1:size(M1,2)];
M2datacol = [max(indexcolumns)+1:size(M2,2)];

%M1datacol = setdiff([1:size(M1,2)],indexcolumns);
%M2datacol = setdiff([1:size(M2,2)],indexcolumns);


matches = rowfind(M1(:,indexcolumns),M2(:,indexcolumns));
out = [];
for i = 1:length(matches)
    if (matches(i) > 0)
        if isempty(outindex)
            out = [out; [M1(i,[indexcolumns M1datacol]) M2(matches(i),M2datacol)]];
        else
            out = [out; [M1(i,[outindex M1datacol]) M2(matches(i),M2datacol)]];
        end
    end
end