function s = tdfread(filename,delimiter,displayopt,readvarnames)
%TDFREAD Read in text and numeric data from tab-delimited file.
%   TDFREAD displays a dialog box for selecting a file, then reads data from
%   the file.  The file should have variable names separated by tabs in the
%   first row, and data values separated by tabs in the remaining rows.
%   TDFREAD creates variables in the workspace, one for each column of the
%   file.  The variable names are taken from the first row of the file.  If a
%   column of the file contains only numeric data in the second and following
%   rows, TDFREAD creates a double variable.  Otherwise, TDFREAD creates a
%   char variable.
%
%   TDFREAD(FILENAME) reads data from the specified file.
%
%   TDFREAD(FILENAME,DELIMITER) uses the specified delimiter in place of tabs.
%   Allowable values of DELIMITER can be any of the following:
%        ' ', '\t', ',', ';', '|' or their corresponding string names 
%        'space', 'tab', 'comma', 'semi', 'bar'
%   'tab' is the default. 
%
%   S = TDFREAD(FILENAME,...) returns a scalar structure S whose fields each
%   contain a variable.
%
%   See also TBLREAD, TEXTSCAN.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.7.4.4 $  $Date: 2006/11/11 22:55:52 $

tab = sprintf('\t');
lf = sprintf('\n');

% set the delimiter
if nargin < 2 || isempty(delimiter)
   delimiter = tab;
else
   switch delimiter
   case {'tab', '\t'}
      delimiter = tab;
   case {'space',' '}
      delimiter = ' ';
   case {'comma', ','}
      delimiter = ',';
   case {'semi', ';'}
      delimiter = ';';
   case {'bar', '|'}
      delimiter = '|';
   otherwise
      delimiter = delimiter(1);
      warning('stats:tdfread:UnrecognizedDelimiter', ...
              ['TDFREAD does not support the character %c as a delimiter.\n' ...
               'This may give bad results.'],delimiter(1));
   end
end

%%% open file
F = 1;
if nargin == 0
   filename = '';
end
if isempty(filename)
   [F,P]=uigetfile('*.*');
   filename = [P,F];
end

if (isequal(F,0)), return, end
fid = fopen(filename,'rt'); % text mode: CRLF -> LF

if fid == -1
   error('stats:tdfread:OpenFailed', ...
         'Unable to open file ''%s''.', filename);
end

% now read in the data
[bigM,count] = fread(fid,Inf);
fclose(fid);
if bigM(count) ~= lf
   bigM = [bigM; lf];
end
bigM = char(bigM(:)');

% remove insignificant whitespace (before and after delimiters or line
% breaks), replace multiple embedded whitespace with a single space, and
% multiple underscores with single underscore.
if delimiter == tab
   matchexpr = {' *\t *' ' *\n+ *' ' +' '_+'};
else
   matchexpr = {['[ \t]*\' delimiter '[ \t]*'] '[ \t]*\n+[ \t]*' '[ \t]+' '_+'};
end
replexpr = {'${delimiter}' '\n' ' ' '_'};
bigM = regexprep(bigM,matchexpr,replexpr);

% find out how many lines are there.
newlines = find(bigM == lf);

% take the first line out from bigM, and put it to line1.
line1 = bigM(1:newlines(1)-1);
if nargin < 4 || readvarnames
   bigM(1:newlines(1)) = [];
   newlines(1) = [];
end

% add a delimiter to the beginning and end of the line
if line1(1) ~= delimiter
   line1 = [delimiter, line1];
end
if line1(end) ~= delimiter
   line1 = [line1, delimiter];
end

% determine varnames
idx = find(line1==delimiter);
nvars = length(idx)-1;
if nargin < 4 || readvarnames
   varnames = cell(1, nvars);
   for k = 1:nvars;
       vn = line1(idx(k)+1:idx(k+1)-1);
       if isempty(vn) % things like ', ,' are already reduced to ',,'
          varnames{k} = strcat('Var',num2str(k,'%d'));
       else
          vn = strrep(vn, ' ', '_');
          varnames{k} = vn;
       end
   end
   varnames = genvarname(varnames);
else
   varnames = strcat({'Var'},num2str((1:nvars)','%d'));
end

nobs = length(newlines);

delimitidx = find(bigM == delimiter);

% check the size validation
if length(delimitidx) ~= nobs*(nvars-1)
   error('stats:tdfread:BadFileFormat',...
         'Requires the same number of delimiters on each line.');
end
if nvars > 1
   delimitidx = (reshape(delimitidx,nvars-1,nobs))';
end

% now we need to re-find the newlines.
newlines = find(bigM(:) == lf);

startlines = [0; newlines(1:nobs-1)];
delimitidx = [startlines, delimitidx, newlines];
fieldlengths = diff(delimitidx,[],2) - 1; fieldlengths = fieldlengths(:);
if any(fieldlengths < 0)
   error('stats:tdfread:BadFileFormat',...
         'Requires the same number of delimiters on each line.');
end
maxlength = max(fieldlengths);
if nargout > 0
    s = struct;
end
for vars = 1:nvars
   x = repmat(' ', nobs,maxlength);
   xempty = false(nobs,1);
   for k = 1:nobs
       v = bigM(delimitidx(k,vars)+1:delimitidx(k,vars+1)-1);
       if isempty(v)
           xempty(k) = true;
       else
           x(k,1:length(v)) = v;
       end
   end
   % Try to convert the strings to create a double variable.  If that
   % succeeds, there may still have been empty strings that were ignored.
   % Treat those as NaNs.  If the conversion failed, there must have been
   % strings that could not be converted.  Leave the variable as char.
   [y,ok] = str2num(x);
   if ok
      if length(y) == size(x,1)
         x = y;
      else
         x = NaN(nobs,1);
         x(~xempty) = y;
      end
   else
      x = deblank(x);
   end
   vname = varnames{vars};
   try
      if nargout == 0
         assignin('base', vname, x);
      else
         s.(vname) = x;
      end
   catch
      warning('stats:tdfread:VarCreateFailed', ...
              'Failed to create variable named %s, skipping column %d from file.', vname, vars);
   end
end

if (nargout == 0) && (nargin<3 || ~isequal(displayopt, 'off'))
   evalin('base', ['whos ', sprintf('%s ',varnames{:})]);
end

