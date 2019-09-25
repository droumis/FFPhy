function display(a)
%@DATASET/DISPLAY Display a CLASSREGTREE object.
%   DISPLAY(T) prints the CLASSREGTREE object T.  DISPLAY is called when a
%   semicolon is not used to terminate a statement.
%
%   See also CLASSREGTREE, CLASSREGTREE/EVAL, CLASSREGTREE/TEST, CLASSREGTREE/PRUNE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:56:17 $

isLoose = strcmp(get(0,'FormatSpacing'),'loose');

objectname = inputname(1);
if isempty(objectname)
   objectname = 'ans';
end

if (isLoose), fprintf('\n'); end
fprintf('%s = \n', objectname);
disp(a)
