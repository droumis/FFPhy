function labels = checklabels(labels,dupFlag)
%CHECKLABELS Validate a list of categorical level text labels.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:24 $

if iscellstr(labels)
    labels = strtrim(labels(:)'); % force cellstr to a column
    if ~all(cellfun(@(s) isvector(s) && (size(s,1) == 1),labels))
        error('stats:categorical:checklabels:InvalidLabels', ...
              '%s must be a character matrix or a cell array of strings.',upper(inputname(1)));
    end
elseif ischar(labels) && (ndims(labels) == 2)
    labels = strtrim(cellstr(labels));
    labels = labels(:)';
else
    error('stats:categorical:checklabels:InvalidLabels', ...
          '%s must be a character matrix or a cell array of strings.',upper(inputname(1)));
end

if any(strcmp(categorical.undefLabel,labels))
    error('stats:categorical:checklabels:UndefinedLabel', ...
          '%s may not contain the string ''%s''.',upper(inputname(1)),categorical.undefLabel);
elseif any(strcmp('',labels))
    error('stats:categorical:checklabels:EmptyLabel', ...
          '%s may not contain the empty string.',upper(inputname(1)));
end

if dupFlag > 0 && length(labels) > 1
    [sortedlevels,ord] = sort(labels);
    d = [true ~strcmp(sortedlevels(2:end),sortedlevels(1:end-1))];
    if ~all(d)
        switch dupFlag
        case 1 % warn if any duplicate labels
            warning('stats:categorical:checklabels:DuplicateLabels', ...
                    'Ignoring duplicate values in %s.',upper(inputname(1)));
            labels = labels(sort(ord(d))); % leave in original order
        case 2 % error if any duplicate labels
            error('stats:categorical:checklabels:DuplicateLabels', ...
                  '%s contains duplicated values.',upper(inputname(1)));
        end
    end
end
