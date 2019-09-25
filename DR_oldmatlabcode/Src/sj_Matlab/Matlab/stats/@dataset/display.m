function display(a)
%DISPLAY Display a dataset array.
%   DISPLAY(DS) prints the dataset array DS, including variable names and
%   observation names (if present).  DISPLAY is called when a semicolon is not
%   used to terminate a statement.
%
%   For numeric or categorical variables that are 2-dimensional and have 3 or
%   fewer columns, DISPLAY prints the actual data.  Otherwise, DISPLAY prints
%   the size and type of each dataset element.
%
%   For character variables that are 2-dimensional and 10 or fewer characters
%   wide, DISPLAY prints quoted strings.  Otherwise, DISPLAY prints the size
%   and type of each dataset element.
%
%   For cell variables that are 2-dimensional and have 3 or fewer columns,
%   DISPLAY prints the contents of each cell (or its size and type if too
%   large). Otherwise, DISPLAY prints the size of each dataset element.
%
%   For time series variables, DISPLAY prints columns for both the time and
%   the data.  If the variable is 2-dimensional and has 3 or fewer columns,
%   DISPLAY prints the actual data.  Otherwise, DISPLAY prints the size and
%   type of each dataset element.
%
%   For other types of variables, DISPLAY prints the size and type of each
%   dataset element.
%
%   See also DATASET, @DATASET/DISP.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:31:35 $

isLoose = strcmp(get(0,'FormatSpacing'),'loose');

objectname = inputname(1);
if isempty(objectname)
   objectname = 'ans';
end

if (isLoose), fprintf('\n'); end
fprintf('%s = \n', objectname);
disp(a)
