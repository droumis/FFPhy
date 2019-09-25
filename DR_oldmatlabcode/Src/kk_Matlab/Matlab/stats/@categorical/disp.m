function disp(a)
%DISP Display a categorical array.
%   DISP(A) prints the categorical array A without printing the array name.
%   In all other ways it's the same as leaving the semicolon off an
%   expression, except that empty arrays don't display.
%
%   See also CATEGORICAL, CATEGORICAL/DISPLAY.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:29 $


isLoose = strcmp(get(0,'FormatSpacing'),'loose');

if (isLoose)
   fprintf('\n');
end
labs = [categorical.undefLabel a.labels];
s = evalc('disp(reshape(labs(a.codes+1),size(a.codes)))');
s(s=='''') = ' ';
disp(s)
