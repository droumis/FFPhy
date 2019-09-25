function newstruct = evaluatefilter2(cellvar, filterString)
%
% out = evaluatefilter(cellvar, filterString)
%
% Return a cell array the same size as the input cell structure.  Each cell
% contains a matrix where the first column is time and the second column is
% 1's and zeros (true or false).
%
% CELLVAR - any cell structure, for example spikes{}{}{}{}.#### There 
%           must be a time field containing a vector of times. The first output column  
%           is time. All other fields
%           are searcheable, and must be the same length as the time field.
% FILTERSTRING - A string containing an expression with the field variables
%               of CELLVAR.  Each field variable must be preceded with a
%               '$'.
%
% Example: index = evaluatefilter2(linbehave,'(($lindist > 10) && ($velocity > 6))')


newstruct = [];
indexstruct = cellfetch(cellvar, '');

[variables, expression] = parsefilterstring(filterString);

if (size(indexstruct.index,2) == 3) %we are search time segments (i.e., trials)
    epochlist = indexstruct.index(:,1:2);
    epochlist = unique(epochlist,'rows');
    for i = 1:size(epochlist,1)
        newstruct{epochlist(i,1)}{epochlist(i,2)} = [-inf 0; inf 0];
    end
end

for i = 1:size(indexstruct.index,1)
    structVar = [];
    eval(['structVar = cellvar', createcellindex(indexstruct.index(i,:)),';']);
    
    for j = 1:length(variables)
        if ( (isfield(structVar,variables{j})) && (isfield(structVar,'time')) && (size(indexstruct.index,2)==2) ) 
            if (size(getfield(structVar,variables{j}),1) ~= size(structVar.time,1))
                error('Filter error: each seach field must be the same length as the time field');
            end
            eval(['tmpresult = ',expression,';']);
            if (size(tmpresult,1) ~= size(structVar.time,1))
                tmpresult = tmpresult';
            end
            eval(['newstruct', createcellindex(indexstruct.index(i,:)),'= [structVar.time tmpresult];']);
            returntype = 1;  
        elseif (isfield(structVar,variables{j}) && isfield(structVar,'timerange') && (size(indexstruct.index,2)==3))
            
            eval(['tmpresult = ',expression,';']);
            if ((tmpresult == 1)||(tmpresult == 0))
                
                newstruct{indexstruct.index(i,1)}{indexstruct.index(i,2)} = [newstruct{indexstruct.index(i,1)}{indexstruct.index(i,2)}; structVar.timerange(1)-.000000000001 0; structVar.timerange(1) tmpresult; structVar.timerange(2) tmpresult; structVar.timerange(2)+.000000000001 0];
            
            elseif length(tmpresult > 1)
                error('Time filter error: for timerange searches, the search parameter must have a single scalar output');
            else
                error('Error in time filter.  Make sure syntax is correct and file is properly searchable.');
            end
              
            returntype = 2;  
        elseif ~isfield(structVar,variables{j})
            error([variables{j},' field does not exist'])
        else
            error('Error in calling timefilter.  Please make sure syntax is correct and file is properly searchable..');
        end
    end
    
end
if (returntype == 2)
    
    for i = 1:size(epochlist,1)
        newstruct{epochlist(i,1)}{epochlist(i,2)} = sortrows(newstruct{epochlist(i,1)}{epochlist(i,2)},1);
    end
end
