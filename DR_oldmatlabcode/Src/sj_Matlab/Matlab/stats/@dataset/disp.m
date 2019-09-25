function disp(a)
%DISP Display a dataset array.
%   DISP(DS) prints the dataset array DS, including variable names and
%   observation names (if present), without printing the dataset name.  In all
%   other ways it's the same as leaving the semicolon off an expression.
%
%   For numeric or categorical variables that are 2-dimensional and have 3 or
%   fewer columns, DISP prints the actual data.  Otherwise, DISP prints the
%   size and type of each dataset element.
%
%   For character variables that are 2-dimensional and 10 or fewer characters
%   wide, DISP prints quoted strings.  Otherwise, DISP prints the size and
%   type of each dataset element.
%
%   For cell variables that are 2-dimensional and have 3 or fewer columns,
%   DISP prints the contents of each cell (or its size and type if too large).
%   Otherwise, DISP prints the size of each dataset element.
%
%   For time series variables, DISP prints columns for both the time and the
%   data.  If the variable is 2-dimensional and has 3 or fewer columns, DISP
%   prints the actual data.  Otherwise, DISP prints the size and type of each
%   dataset element.
%
%   For other types of variables, DISP prints the size and type of each
%   dataset element.
%
%   See also DATASET, @DATASET/DISPLAY.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:31:34 $

isLoose = strcmp(get(0,'FormatSpacing'),'loose');
maxWidth = get(0,'CommandWindowSize'); maxWidth = maxWidth(1);

if (isLoose), fprintf('\n'); end

if (a.nobs > 0) && (a.nvars > 0)
    dsPad = repmat(' ', a.nobs+1, 4);
    varPad = repmat(' ', a.nobs, 2);
    if isempty(a.obsnames)
        dsChars = char(zeros(a.nobs+1, 0));
    else
        dsChars = [dsPad strvcat(' ',char(a.obsnames))];
    end
    for i = 1:a.nvars
        name = a.varnames{i};
        var = a.data{i};
        
        if ischar(var)
            if (ndims(var) == 2) && (size(var,2) <= 10)
                % Display individual strings for a char variable that is 2D and no
                % more than 10 chars.
                varChars = var;
            else
                % Otherwise, display a description of the chars.
                sz = size(var);
                if ndims(var) == 2
                    szStr = ['[1' sprintf('x%d',sz(2:end))];
                else
                    szStr = ['[' sprintf('%d',sz(2)) sprintf('x%d',sz(3:end))];
                end
                varChars = repmat([szStr ' char]'],sz(1),1);
            end
            
        % Display both the time and the data for a time series variable
        elseif isa(var,'timeseries')
            if (ndims(var.data) == 2) && (size(var.data,2) <= 3)
                % Display the individual data if the var is 2D and no more than 3 columns.
                varChars = [num2str(var.time) repmat('  ',a.nobs,1) num2str(var.data)];
            else
                % Otherwise display a description of the data
                sz = size(var.data);
                if ndims(var.data) == 2
                    szStr = ['[1' sprintf('x%d',sz(2:end))];
                else
                    szStr = ['[' sprintf('%d',sz(2)) sprintf('x%d',sz(3:end))];
                end
                varChars = [num2str(var.time) repmat('  ',a.nobs,1) repmat([szStr ' ' class(var.data) ']'],sz(1),1)];
            end
                
        else
            % Display the individual data if the var is 2D and no more than 3 columns.
            if ~isempty(var) && (ndims(var) == 2) && (size(var,2) <= 3)
                if isnumeric(var)
                    varChars = num2str(var);
                elseif islogical(var)
                    % Display the logical values using meaningful names.
                    strs = ['false'; 'true '];
                    w1 = size(strs,2); w2 = size(varPad,2);
                    varChars = repmat(' ',size(var,1),(size(var,2)-1)*(w1+w2));
                    for j = 1:size(var,2)
                        varChars(:,(j-1)*(w1+w2)+(1:w1)) = strs(var(:,j)+1,:);
                    end
                elseif isa(var,'categorical')
                    % Build the output one column at a time, since the char method reshapes
                    % to a single column.
                    varChars = char(zeros(a.nobs,0));
                    for j = 1:size(var,2)
                        if j > 1, varChars = [varChars varPad]; end
                        varChars = [varChars char(var(:,j))];
                    end
                elseif iscell(var)
                    % Display the contents of each cell, or a description, as
                    % the built-in cell display method sees fit.
                    varStr = evalc('disp(var)');
                    loc = [0 find(varStr==10)];
                    len = diff(loc);
                    varChars = repmat(' ',size(var,1),max(len)-1);
                    for j = 1:size(var,1)
                        celChars = strtrim(varStr(loc(j)+1:loc(j+1)-1));
                        varChars(j,1:length(celChars)) = celChars;
                    end
                else
                    % Display a description of each dataset element.
                    sz = size(var);
                    szStr = ['[1' sprintf('x%d',sz(2:end))];
                    varChars = repmat([szStr ' ' class(var) ']'],sz(1),1);
                end

            % Either the variable is not 2D, or it's empty, or it's too wide
            % to show. Display a description of each dataset element.
            else
                sz = size(var);
                if ndims(var) == 2
                    szStr = ['[1' sprintf('x%d',sz(2:end))];
                else
                    szStr = ['[' sprintf('%d',sz(2)) sprintf('x%d',sz(3:end))];
                end
                varChars = repmat([szStr ' ' class(var) ']'],sz(1),1);
            end
        end
        varChars = strvcat(name,varChars);
        
        % If we've reached the right margin, display the output built up so far, and then
        % restart for display starting at the left margin.
        if size(dsChars,2) + size(dsPad,2) + size(varChars,2) > maxWidth
            disp(dsChars)
            fprintf('\n');
            if (isLoose)
                fprintf('\n');
            end
            if isempty(a.obsnames)
                dsChars = char(zeros(a.nobs+1, 0));
            else
                dsChars = [dsPad strvcat(' ',char(a.obsnames))];
            end
        end
        dsChars = [dsChars dsPad varChars];
    end
else
    dsChars = sprintf('[empty %d-by-%d dataset]',a.nobs,a.nvars);
end
disp(dsChars)
if (isLoose), fprintf('\n'); end
