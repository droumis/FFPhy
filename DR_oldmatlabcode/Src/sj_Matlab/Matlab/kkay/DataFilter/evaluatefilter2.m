function newstruct = evaluatefilter2(cellvar, filterString)
%
% out = evaluatefilter(cellvar, filterString)
%
% Return a cell array the same size as the input cell structure.  Each cell
% contains a matrix where the first column is time and the second column is
% 1's and zeros (true or false).
%
% CELLVAR - any cell structure, for example spikes{}{}{}{}.#### There must
%           be a time field containing a vector of times.  All other fields
%           are searcheable, and must be the same length as the time field.
% FILTERSTRING - A string containing an expression with the field variables
%               of CELLVAR.  Each field variable must be preceded with a
%               '$'.
%
% Example: index = evaluatefilter2(linbehave,'(($lindist > 10) && ($velocity > 6))')


newstruct = [];
indexstruct = cellfetch(cellvar, '');
[variables, expression] = parsefilterstring(filterString);
for i = 1:size(indexstruct.index,1) 
    structVar = [];
    eval(['structVar = cellvar', createcellindex(indexstruct.index(i,:)),';']);
    for j = 1:length(variables)
        if isfield(structVar,variables{j})
            if (size(getfield(structVar,variables{j}),1) ~= size(structVar.time,1))
                error('Filter error: each seach field must be the same length as the time field');
            end
        else
            error([variables{j},' field does not exist'])
        end
    end
    eval(['tmpresult = ',expression,';']);
    % if the data are column vectors, transpose them
    if (size(tmpresult,2) > size(tmpresult,1))
	tmpresult = tmpresult';
    end
    if (size(structVar.time,2) > size(structVar.time,1))
	structVar.time = structVar.time';
    end
    eval(['newstruct', createcellindex(indexstruct.index(i,:)),'= [structVar.time tmpresult];']);
end
end

