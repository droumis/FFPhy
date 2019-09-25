classdef dataset
%DATASET Create a dataset array.
%   DS = DATASET(VAR1, VAR2, ...) creates a dataset array DS from the
%   workspace variables VAR1, VAR2, ... .  All variables must have the same
%   number of rows.
%
%   DS = DATASET(..., {VAR,'name'}, ...) creates a dataset variable named
%   'name' in DS.  Dataset variable names must be valid MATLAB identifiers,
%   and unique.
%
%   DS = DATASET(..., {VAR,'name1',...,'name_M'}, ...), where VAR is an
%   N-by-M-by-P-by-... array, creates M dataset variables in DS, each of size
%   N-by-P-by-..., with names 'name1', ..., 'name_M'.
%
%   DS = DATASET(..., 'varnames', {'name1', ..., 'name_M'}) creates dataset
%   variables that have the specified variable names.  The names must be valid
%   MATLAB identifiers, and unique.  You may not provide both the 'varnames'
%   parameter and names for individual variables.
%
%   DS = DATASET(..., 'obsnames', {'name1', ..., 'name_N'}) creates a dataset
%   array that has the specified observation names.  The names need not be
%   valid MATLAB identifiers, but must be unique.
%
%   Dataset arrays can contain variables that are built-in types, or objects that
%   are arrays and support standard MATLAB parenthesis indexing of the form
%   var(i,...), where i is a numeric or logical vector that corresponds to
%   rows of the variable.  In addition, the array must implement a SIZE method
%   with a DIM argument, and a VERTCAT method.
%
%   DS = DATASET('File',FILENAME, ...) creates a dataset array from
%   column-oriented data in a text file.  Variable names are taken from the
%   first row of the file. You may also specify a delimiter character using
%   the 'Delimiter' parameter name/value pair.  This may be any delimiter
%   accepted by the TDFREAD function.
%
%   DS = DATASET('File',FILENAME,'Format',FORMAT, ...)  creates a dataset
%   array from column-oriented data in a text file, where FORMAT is a format
%   string as accepted by the TEXTSCAN function.  You may also specify any of
%   the parameter name/value pairs accepted by the TEXTSCAN function.
%   Variable names are taken from the first row of the file.
%
%   DS = DATASET('XLSFile',XLSFILENAME, ...) creates a dataset array from
%   column-oriented data in an Excel spreadsheet file.  You may also specify
%   the 'Sheet' and 'Range' parameter name/value pairs as accepted by the
%   XLSREAD function.  Variable names are taken from the first row of the
%   spreadsheet.
%
%   When reading from a text or spreadsheet file, the 'ReadVarNames' parameter
%   name/value pair determines whether or not the first row of the file is
%   treated as variable names.  Specify as a logical value (default true).
%
%   When reading from a text or spreadsheet file, the 'ReadObsNames' parameter
%   name/value pair determines whether or not the first column of the file is
%   treated as observation names.  Specify as a logical value (default false).
%   If the 'ReadVarNames' and 'ReadObsNames' parameter values are both true,
%   the name in the first column of the first row of the file is saved as the
%   first dimension name for the dataset.
%
%   Reading from a text or spreadsheet file creates scalar-valued dataset
%   variables, i.e., one variable from each column in the file.  The dataset
%   variables that are created are either double-valued, if the entire column
%   is numeric, or string-valued, i.e. a cell array of strings, if any element
%   in a column is not numeric.
%
%   See also DATASET/SET, DATASET/GET, GENVARNAME, TDFREAD, TEXTSCAN, XLSREAD.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:31:32 $

    properties(GetAccess='private', SetAccess='private')
        ndims = 2;
        nobs = 0;
        obsnames = {};
        nvars = 0;
        varnames = cell(1,0); % these can never be "truly" empty
        data = {};
        
        % 'Properties' will also appear to contain 'VarNames' and 'ObsNames'
        props = struct('Description',{''}, ...
                       'Units',{{}}, ...
                       'DimNames',{{'Observations' 'Variables'}}, ...
                       'UserData',{[]});
    end
    
    methods
        function a = dataset(varargin)
        
        varnames = repmat({''},1,nargin);
        hadExplicitNames = false;
        argCnt = 0;
        
        if nargin == 0
            % nothing to do
            return
        end
        
        % Set these as a first guess, they may grow if some vars are
        % created by splitting up arrays, or shrink if there are
        % name/value pairs.
        a.nvars = nargin;
        a.data = cell(1,nargin);

        varCnt = 0;
        while argCnt < nargin
            argCnt = argCnt + 1;
            arg = varargin{argCnt};
            if isstring(arg)
                % start of param name/value pairs
                break
            elseif iscell(arg) && ~isscalar(arg) && isvector(arg) && (size(arg,1)==1)
                if (numel(arg)==2) && isstring(arg{2})
                    % {var,name}
                    varCnt = varCnt + 1;
                    a.data{varCnt} = arg{1};
                    varnames{varCnt} = arg{2};
                    hadExplicitNames = true;
                elseif (numel(arg)>2) && all(cellfun(@isstring,arg(2:end)))
                    % {var,name1,name2,...}
                    var = arg{1};                    
                    names = arg(2:end);
                    if length(names) ~= size(var,2)
                        error('stats:dataset:dataset:IncorrectNumberVarnames', ...
                              'Must have one variable name for each column when creating multiple variables from an array.');
                    end
                    szOut = size(var);
                    szOut(2) = []; if isscalar(szOut), szOut(2) = 1; end
                    ncols = size(var,2);
                    a.nvars = a.nvars + ncols-1;
                    varnames = [varnames repmat({''},1,ncols-1)];
                    a.data = [a.data cell(1,ncols-1)];
                    for j = 1:ncols
                        varCnt = varCnt + 1;
                        a.data{varCnt} = reshape(var(:,j,:),szOut);
                        varnames{varCnt} = names{j};
                    end
                    hadExplicitNames = true;
                else
                    % false alarm -- cell-valued var without name
                    varCnt = varCnt + 1;
                    a.data{varCnt} = arg;
                    name = inputname(argCnt);
                    if isempty(name)
                        name = strcat('Var',num2str(varCnt,'%d'));
                    end
                    varnames{varCnt} = name;
                end
            elseif isstruct(arg) && isscalar(arg)
                if any(diff(structfun(@(x)size(x,1),arg)))
                    error('stats:dataset:dataset:UnequalFieldLengths', ...
                          'Fields in a scalar structure must have the same number of rows.');
                end
                names = fieldnames(arg);
                ncols = length(names);
                a.nvars = a.nvars + ncols-1;
                varnames = [varnames repmat({''},1,ncols-1)];
                a.data = [a.data cell(1,ncols-1)];
                for j = 1:ncols
                    varCnt = varCnt + 1;
                    a.data{varCnt} = arg.(names{j});
                    varnames{varCnt} = names{j};
                end
                hadExplicitNames = true;
            else
                % var without name
                varCnt = varCnt + 1;
                a.data{varCnt} = arg;
                name = inputname(argCnt);
                if isempty(name)
                    name = strcat('Var',num2str(varCnt,'%d'));
                end
                varnames{varCnt} = name;
            end
            try
                nrows = size(a.data{varCnt},1);
            catch
                error('stats:dataset:dataset:SizeMethodFailed', ...
                      'SIZE(VAR,1) for variable %d failed with the following error:\n\n%s.',varCnt,lasterr);
            end
            if argCnt == 1
                a.nobs = nrows;
            elseif ~isequal(nrows,a.nobs)
                error('stats:dataset:dataset:UnequalVarLengths', ...
                      'All variables must have the same number of rows.');
            end

        end % while argCnt < nargin, processing individual vars
        
        a.nvars = varCnt;
        a.data = a.data(1:varCnt);
        varnames = varnames(1:varCnt);
            
        if argCnt < nargin
            pnames = {'file' 'xlsfile' 'varnames'  'obsnames'};
            dflts =  {    []        []         []          []};
            [eid,errmsg,fileArg,xlsfileArg,varnamesArg,obsnamesArg,otherArgs] ...
                      = dfswitchyard('statgetargs', pnames, dflts, varargin{argCnt:end});
            if ~isempty(eid)
                error(sprintf('stats:dataset:dataset:%s',eid),errmsg);
            end
            
            if ~isempty(fileArg)
                if argCnt > 1
                    error('stats:dataset:dataset:FileAndData', ...
                          'You cannot specify the ''file'' parameter and individual variables.');
                end
                a = readFile(a,fileArg,otherArgs);
            elseif ~isempty(xlsfileArg)
                if argCnt > 1
                    error('stats:dataset:dataset:XLSFileAndData', ...
                          'You cannot specify the ''xlsfile'' parameter and individual variables.');
                end
                a = readXLSFile(a,xlsfileArg,otherArgs);
            else
                if ~isempty(otherArgs)
                    error('stats:dataset:dataset:UnrecognizedParams', ...
                          'The parameter %s is unrecognized or not legal in this context.',otherArgs{1});
                end
            end
            
            if ~isempty(varnamesArg)
                if hadExplicitNames
                    error('stats:dataset:dataset:VarNamesAndVarNamesParam', ...
                          'You cannot specify the ''VarNames'' parameter and individual variable names.');
                else
                    varnames = varnamesArg;
                end
            end
            if ~isempty(obsnamesArg)
                a = setobsnames(a,obsnamesArg);
            end
        end % if argCnt < nargin,processing name/value pairs
        
        % Varnames may be empty because we had no vars, or becasue we read
        % from a file.  In either case, no need to set them.
        if ~isempty(varnames)
            a = setvarnames(a,varnames);
        end
        
        end % dataset constructor
    end % methods block
    
    methods(Visible = false)
        function b = fieldnames(a)
            b = a.varnames(:);
            b{end+1} = 'Properties';
        end
        
        % Methods that we inherit from opaque, but do not want
        function a = fields(varargin),          throwUndefinedError; end
        function a = toChar(varargin),          throwUndefinedError; end

        function a = permute(varargin),         throwUndefinedError; end
        function a = transpose(varargin),       throwUndefinedError; end
        function a = ctranspose(varargin),      throwUndefinedError; end
        function a = reshape(varargin),         throwUndefinedError; end
        function a = tril(varargin),            throwUndefinedError; end
        function a = triu(varargin),            throwUndefinedError; end
        function a = diag(varargin),            throwUndefinedError; end

        function [a,b] = sort(varargin),        throwUndefinedError; end
        function [a,b] = ismember(varargin),    throwUndefinedError; end
        function [a,b] = setdiff(varargin),     throwUndefinedError; end
        function [a,b,c] = setxor(varargin),    throwUndefinedError; end
        function [a,b,c] = intersect(varargin), throwUndefinedError; end
        
    end % methods block
end % classdef


function tf = isstring(s) % require a row of chars, or possibly ''
tf = ischar(s) && ( (isvector(s) && (size(s,1) == 1)) || all(size(s)==0) );
end

function throwUndefinedError
st = dbstack;
name = strread(st(2).name,'dataset.%s');
error('stats:dataset:UndefinedFunction', ...
      'Undefined function or method ''%s'' for input arguments of type ''dataset''.',name{1});
end



%-----------------------------------------------------------------------------
function a = readFile(a,file,args)

pnames = {'readvarnames' 'readobsnames' 'delimiter' 'format' 'headerlines'};
dflts =  {          true          false          []       []            []};
[eid,errmsg,readvarnames,readobsnames,delimiter,format,headerlines,otherArgs] ...
                   = dfswitchyard('statgetargs', pnames, dflts, args{:});
if ~isempty(eid)
    error(sprintf('stats:dataset:dataset:%s',eid),errmsg);
end
readobsnames = onOff2Logical(readobsnames,'ReadObsNames');
readvarnames = onOff2Logical(readvarnames,'ReadVarNames');

if ~isempty(format)
    if isempty(delimiter), delimiter = ' '; end
    if isempty(headerlines), headerlines = 0; end
    fid = fopen(file);
    if fid == -1
        error('stats:dataset:dataset:OpenFailed', ...
              'Unable to open the file %s for reading.',file);
    end
    if readvarnames
        % Search for any of the allowable conversions: a '%', followed
        % optionally by '*', followed optionally by 'nnn' or by 'nnn.nnn',
        % followed by one of the type specifiers or a character list or a
        % negative character list.  Keep the '%' and the '*' if it's there,
        % but replace everything else with the 'q' type specifier.
        specifiers = '(n|d8|d16|d32|d64|d|u8|u16|u32|u64|u|f32|f64|f|s|q|c|\[\^?[^\[\%]*\])';
        vnformat = regexprep(format,['\%([*]?)([0-9]+(.[0-9]+)?)?' specifiers],'%$1q');
        varnames = textscan(fid,vnformat,1,'delimiter',delimiter,'headerlines',headerlines,otherArgs{:});
        % If a textscan was unable to read a varname, the corresponding cell
        % contains an empty cell.  Remove those.  This happens when trying to
        % put delimiters in the format string, because they're just read as
        % part of the string.
        varnames = varnames(~cellfun('isempty',varnames));
        % Each cell in varnames contains another 1x1 cell containing a string,
        % get those out.
        varnames = cellfun(@(c) c{1}, varnames,'UniformOutput',false);
        headerlines = 0; % just skipped them
    end
    raw = textscan(fid,format,'delimiter',delimiter,'headerlines',headerlines,otherArgs{:});
    fclose(fid);
    if readvarnames
        if numel(varnames) ~= numel(raw)
            error('stats:dataset:dataset:ReadVarnamesFailed', ...
                  ['The number of variable names read from %s does not match the number ' ...
                   'of\ndata columns read.  You may have specified the format string, ' ...
                   'delimiter,\nor the number of header lines incorrectly.'],file);
        end
        varnames = genuniquenames(genvalidnames(varnames),1);
    else
        varnames = strcat({'Var'},num2str((1:numel(raw))','%d'));
    end
    s = cell2struct(raw(:),varnames,1);
else
    if ~isempty(headerlines) || ~isempty(otherArgs)
        if ~isempty(headerlines), arg = 'HeaderLines'; else, arg = otherArgs{1}; end
        error('stats:dataset:dataset:UnrecognizedParams', ...
              'The ''%s'' parameter is unrecognized or not legal in this context.',arg);
    end
    s = tdfread(file,delimiter,false,readvarnames);
    varnames = fieldnames(s);
    for i = 1:length(varnames)
        name = varnames{i};
        if ischar(s.(name))
            s.(varnames{i}) = cellstr(s.(name));
        end
    end
end

varlen = unique(structfun(@(x)size(x,1),s));
if ~isscalar(varlen)
    errMsg = 'Variable lengths must all be the same.';
    if ~isempty(format)
        errMsg = [errMsg '  You may have specified the\nformat string, delimiter, ' ...
                  'or number of header lines incorrectly.'];
    end
    error('stats:dataset:dataset:UnequalVarLengthsFromFile', errMsg);
end
a.nobs = varlen;
if readobsnames
    obsnames = cellstr(s.(varnames{1}));
    s = rmfield(s,varnames{1});
    dimnames = a.props.DimNames;
    dimnames{1} = varnames{1};
    varnames(1) = [];
    a = setobsnames(a,obsnames);
    a = setdimnames(a,dimnames);
end
a.nvars = length(varnames);
a.data = struct2cell(s); a.data = a.data(:)';
a = setvarnames(a,varnames(:)');

end % readFile function


%-----------------------------------------------------------------------------
function a = readXLSFile(a,xlsfile,args)

pnames = {'readvarnames' 'readobsnames' 'sheet' 'range'};
dflts =  {          true          false      ''      ''};
[eid,errmsg,readvarnames,readobsnames,sheet,range] ...
                   = dfswitchyard('statgetargs', pnames, dflts, args{:});
if ~isempty(eid)
    error(sprintf('stats:dataset:dataset:%s',eid),errmsg);
end
readobsnames = onOff2Logical(readobsnames,'ReadObsNames');
readvarnames = onOff2Logical(readvarnames,'ReadVarNames');

[numeric,txt,raw] = xlsread(xlsfile,sheet,range);
if isempty(numeric) && isempty(txt)
    return
end
clear numeric txt

if readvarnames
    varnames = raw(1,:);
    if ~iscellstr(varnames)
        varnames = cellfun(@convert2Str, varnames, 'UniformOutput',false);
    end
    raw(1,:) = [];
else
    varnames = strcat({'Var'},num2str((1:size(raw,2))','%d'));
end

a.nobs = size(raw,1);
if readobsnames
    obsnames = raw(:,1);
    if ~iscellstr(obsnames)
        obsnames = cellfun(@convert2Str, obsnames, 'UniformOutput',false);
    end
    dimnames = a.props.DimNames;
    dimnames{1} = varnames{1};
    varnames(1) = [];
    raw(:,1) = [];
    a = setobsnames(a,obsnames);
    a = setdimnames(a,dimnames);
end

a.nvars = length(varnames);
a.data = cell(1,a.nvars);
a = setvarnames(a,varnames(:)');

for j = 1:a.nvars
    if all(cellfun(@isnumeric,raw(:,j)))
        a.data{j} = cell2mat(raw(:,j));
    else
        a.data{j} = cellfun(@convert2Str, raw(:,j), 'UniformOutput',false);
    end
end

end % readXLSFile function


%-----------------------------------------------------------------------------
function s = convert2Str(n)
if isnumeric(n)
    if isnan(n)
        % A numeric NaN means the cell was empty, make that the empty string
        s = '';
    else
        s = num2str(n);
    end
elseif ischar(n)
    s = n;
else
    error('stats:dataset:dataset:UnexpectedClass', ...
          'Unable to convert value of class %s to char.',class(n));
end
end % function
    

%-----------------------------------------------------------------------------
function arg = onOff2Logical(arg,str)

if ischar(arg)
    if strcmpi(arg,'on')
        arg = true;
    elseif strcmpi(arg,'off')
        arg = false;
    else
        error(['stats:dataset:dataset:Invalid' str], ...
              'The ''%s'' parameter must be ''on'', ''off'', or a logical scalar.',str);
    end
elseif islogical(arg)
    % leave it alone
else
    error(['stats:dataset:dataset:Invalid' str], ...
           'The ''%s'' parameter must be ''on'', ''off'', or a logical scalar.',str);
end
end % function
