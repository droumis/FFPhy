function display(a)
%DISPLAY Display a categorical array.
%   DISP(A) prints the categorical array A.  DISPLAY is called when a
%   semicolon is not used to terminate a statement.
%
%   See also CATEGORICAL, CATEGORICAL/DISP.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:30 $

isLoose = strcmp(get(0,'FormatSpacing'),'loose');

objectname = inputname(1);
if isempty(objectname)
    objectname = 'ans';
end

if (isLoose)
    fprintf('\n');
end
if isempty(a)
    sz = size(a);
    if ndims(a) == 2
        fprintf('Empty %s matrix: %d-by-%d\n',class(a),sz(1),sz(2));
    else
        fprintf('Empty %s array: %d',class(a),sz(1));
        fprintf('-by-%d',sz(2:end));
        fprintf('\n');
    end
elseif ndims(a) == 2
    fprintf('%s = \n', objectname);
    disp(a)
else
    if (isLoose)
        fprintf('\n');
    end
    labs = [categorical.undefLabel a.labels];
    % Let the builtin cellstr display method do the real work, then remove
    % quotes, and look for the page headers, things like '(:,:,1) =', and
    % replace them with 'objectname(:,:,1) ='
    s = evalc('disp(reshape(labs(a.codes+1),size(a.codes)))');
    s(s=='''') = ' ';
    s = regexprep(s,'\([0-9:,]+\) =', [objectname '$0']);
    disp(s)
end
