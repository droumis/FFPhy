function s = summary(a)
%SUMMARY Print summary of a dataset array.
%   SUMMARY(A) prints a summary of a dataset array and the variables that it
%   contains.
%
%   S = SUMMARY(A) returns a scalar structure S that contains a summary of the
%   dataset A and the variables that A contains. S contains the following
%   fields:
%
%    Description   A character array containing the dataset description.
%    Variables     A structure array with one element for each dataset
%                  variable in A.  Each element has the following fields:
%       Name       A character string containing the name of the variable.
%       Units      A character string containing the variable's units.
%       Size       A numeric vector containing the size of the variable.
%       Class      A character string containing the class of the variable.
%       Data       A scalar structure containing the following fields:
%          (for numeric variables)
%             Probabilities   A numeric vector containing the probabilities 
%                             [0.0 .25 .50 .75 1.0] and NaN (if any are present
%                             in the corresponding dataset variable).    
%             Quantiles       A numeric vector containing the values that
%                             correspond to 'Probabilities' for the
%                             corresponding dataset variable, and a count of
%                             NaNs (if any are present).
%          (for logical variables)
%             Values          The logical vector [true false].
%             Counts          A numeric vector of counts for each logical value.
%          (for categorical variables)
%             Levels          A cell array containing the labels for each level
%                             of the corresponding dataset variable.
%             Counts          A numeric vector of counts for each level.
%
%   'Data' is empty if variable is not numeric, categorical, or logical.  If a
%   dataset variable has more than one column, then the corresponding
%   'Quantiles' or 'Counts' field is a matrix or an array.
%
%   See also DATASET/GET, DATASET/SET.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:31:51 $

descr = get(a,'Description');
units = get(a,'Units');
varnames = get(a,'VarNames');
    
if nargout < 1
    isLoose = strcmp(get(0,'FormatSpacing'),'loose');

    if (isLoose), fprintf('\n'); end
    
    if ~isempty(descr)
        disp(descr); fprintf('\n');
        if (isLoose), fprintf('\n'); end
    end

    for j = 1:a.nvars
        var_j = a.data{j};
        sz_j = size(var_j);
        szStr = [sprintf('%d',sz_j(1)) sprintf('x%d',sz_j(2:end))];
        if ~isempty(units) && ~isempty(units{j})
            unitsStr = [', Units = ' units{j}];
        else
            unitsStr = '';
        end
        if iscellstr(var_j)
            disp(sprintf('%s: [%s cell string%s]',varnames{j},szStr,unitsStr));
        else
            disp(sprintf('%s: [%s %s%s]',varnames{j},szStr,class(var_j),unitsStr));
        end
        if (sz_j(1) > 0) && (length(sz_j) == 2)
            if isfloat(var_j)
                [q,p,labs] = ctsQuantiles(var_j);
            elseif isinteger(var_j)
                [q,p,labs] = intQuantiles(var_j);
            elseif islogical(var_j)
                [q,tf,labs] = logicalCounts(var_j);
            elseif isa(var_j,'categorical')
                [q,labs] = summary(var_j);
            else
                q = []; labs = {};
            end
            c = [labs num2cell(q)];
            if size(var_j,2) == 1, c = c'; end
            str = evalc('disp(c)');
        else
            str = ' ';
        end
        str = str(1:end-1); % remove trailing newline
        % First find brackets containing numbers in any format, and preceeded by
        % whitespace -- those are the counts/quantiles.  Replace those enclosing
        % brackets with spaces.  Then replace all quotes with spaces.
        str = regexprep(str,'(\s)\[([^\]]+)\]','$1 $2 ');
        str = regexprep(str,'''',' ');
        disp(str);
    end
    
else
    data = cell(size(a.data));
    for j = 1:a.nvars
        var_j = a.data{j};
        if isfloat(var_j)
            [q,p] = ctsQuantiles(var_j);
            data{j} = struct('Probabilities',p, 'Quantiles',q);
        elseif isinteger(var_j)
            [q,p] = intQuantiles(var_j);
            data{j} = struct('Probabilities',p, 'Quantiles',q);
        elseif islogical(var_j)
            [q,tf] = logicalCounts(var_j);
            data{j} = struct('Values',tf, 'Counts',q);
        elseif isa(var_j,'categorical')
            [q,labs] = summary(var_j);
            data{j} = struct('Levels',{labs}, 'Counts',q);
        end
    end
    t = struct('Name', varnames, ...
               'Units', units, ...
               'Size', cellfun(@size, a.data, 'UniformOutput',false), ...
               'Class', cellfun(@class, a.data, 'UniformOutput',false), ...
               'Data', data);
    s = struct('Description',descr,'Variables',t);
end


%-----------------------------------------------------------------------------
function [q,p,labs] = ctsQuantiles(x)
labs = {'min'; '1st Q'; 'median'; '3rd Q'; 'max'};
p = [0; .25; .5; .75; 1];
q = quantile(x,p,1);
nnans = sum(isnan(x),1);
if any(nnans(:) > 0)
    labs = [labs; 'NaNs'];
    p = [p; NaN];
    q = [q; nnans];
end


%-----------------------------------------------------------------------------
function [p,q,labs] = intQuantiles(x)
labs = {'min'; '1st Q'; 'median'; '3rd Q'; 'max'};
p = [0; .25; .5; .75; 1];
x = sort(x,1);
j0 = [.25 .5 .75]*size(x,1) + .5; j1 = floor(j0); j2 = ceil(j0);
a = x(j1,:); b = x(j2,:);
quartiles = a;
w = j0 - j1;
for i = find(w)
    quartiles(i,:) = a(i,:) + w(i)*(b(i,:)-a(i,:));
    k = (sign(a(i,:)) ~= sign(b(i,:)));
    quartiles(i,k) = (1-w(i))*a(i,k) + w(i)*b(i,k);
end
q = [x(1,:); quartiles; x(n,:)];
szOut = size(x); szOut(1) = size(q,1);
q = reshape(q,szOut);


%-----------------------------------------------------------------------------
function [q,tf,labs] = logicalCounts(x)
labs = {'true'; 'false'};
tf = [true; false];
q = [sum(x,1); sum(1-x,1)];

