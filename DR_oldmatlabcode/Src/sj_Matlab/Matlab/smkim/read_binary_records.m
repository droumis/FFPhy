function [recstruct, count, header] = read_binary_records(filename,eoh,recformat,nrecs,varargin)
%READ_BINARY_RECORDS Read C-style struct records from a file.  
%   RECSTRUCT = READ_BINARY_RECORDS(FILENAME,EOH,RECFORMAT,NRECS) reads packed
%   binary data that is formatted in C-style struct records and returns a MATLAB
%   struct. THIS FUNCTION ONLY WORKS ON PLATFORMS THAT HAVE 8 BITS PER BYTE!!!
%
%   EOH is a regular expression for the character sequence that marks the end of
%   the file header; set this to the empty string '' if the file has no header.
%   See the documentation for REGEXP for details on how to construct a regular
%   expression that matches what you intend. It is very important to include
%   any terminating newline character (\n) if there is one at the end of the
%   header.
%
%   RECFORMAT is a struct array that specifies the order and format in which the
%   data are packed. It must have the following fields: 
%      RECFORMAT(i).name: (string) name of the ith field of the record.  
%      RECFORMAT(i).type: string code that specifies the datatype in this field;
%                         field; allowed values are 'schar', 'uchar', 'int8',
%                         'int16', 'int32', 'int64', 'float32', 'float64'. Note
%                         that this function uses hard-coded values for the
%                         sizes of these types that are only correct if 1 byte
%                         = 8 bits. This function will fail on platforms that
%                         do not use octets for bytes!
%     RECFORMAT(i).count: number of elements of this datatype that appear in
%                         the ith field.
%  
%   NRECS is the number of records to read from the file, starting from the
%   beginning. READ_BINARY_RECORDS will read as many complete records as it can
%   find; supplying a value Inf will read to the end of the file.
%
%   RECSTRUCT = READ_BINARY_RECORDS(FILENAME,EOH,RECFORMAT,NUMRECS,ENDDIANNESS)
%   reads the data according to the byte order specified by the string code
%   ENDIANNESS; see documentation for FOPEN for allowed values. By default, this
%   function assumes that the file has the same enddianness as the system on
%   which Matlab is being run.
%
%   [RECSTRUCT, COUNT] = READ_BINARY_RECORDS(...) returns an integer COUNT,
%   which gives the number of records that were successfully read (less than or
%   equal to NRECS).
%
%   [RECSTRUCT, COUNT, HEADER] = READ_BINARY_RECORDS(...) returns an additional
%   char array HEADER, containing the characters at the beginning of the file
%   leading up to and including the EOH character sequence. 
%
%   Example: Suppose that you have a binary file named 'myfile.dat' in which the
%   header is terminated by the character sequence "%%ENDHEADER\n". This file
%   contains records of the following form: a uint32, which you want to call a
%   'timestamp', and 160 uint16 values, which you want to call 'samples'. To
%   read this data into a struct, you would do the following:
%     myformat(1).name = 'timestamp';
%     myformat(1).type = 'uint32';
%     myformat(1).count = 1;
%     myformat(2).name = 'samples';
%     myformat(2).type = 'uint16';
%     myformat(2).count = 160;
%     read_binary_records('myfile.dat','%%ENDHEADER\n',myformat,Inf)
%     
% Written by smk, 2009 February 06.
%

% Note: These hard-coded values are only valid for systems which have 8 bits per
% byte. Proceed with caution. For ease of implementation and code
% maintainability, only a subset of the formats supported by FREAD are supported
% here.  
ALLOWED_TYPES = cell2struct({ ...
    'schar'   , 1; ...
    'uchar'   , 1; ...
    'int8'    , 1; ...
    'int16'   , 2; ...
    'int32'   , 4; ...
    'int64'   , 8; ...
    'uint8'   , 1; ...
    'uint16'  , 2; ...
    'uint32'  , 4; ...
    'uint64'  , 8; ...
    'float32' , 4; ...
    'float64' , 8 }, ...
    {'type', 'nbytes'}, 2);
% helper function associated with this struct
SIZE_OF = @(type_str) ...
    ALLOWED_TYPES(find(strcmp(type_str,{ALLOWED_TYPES(:).type}))).nbytes;
REQUIRED_RECFORMAT_FIELDS = {'count'; 'name'; 'type'};

if (length(varargin) > 1)
  error('too many arguments');
elseif (length(varargin) == 1)
  if ischar(varargin{1})
    machineformat = varargin{1};
  elsea
    error('endianness must be a string code; see help for FOPEN');
  end
else
  machineformat = 'native';
end

if ~ischar(eoh)
  error('eoh must be a string (end-of-header character sequence)');
end

if ~isstruct(recformat) || isempty(recformat)
  error('recformat must be a non-empty struct');
end
if (length(fieldnames(recformat)) ~= 3) || ...
    ~all(isfield(recformat,REQUIRED_RECFORMAT_FIELDS))
  error('recformat has incorrect fields');
end
for i = 1:length(recformat)
  if ~ischar(recformat(i).type) || ...
      ~any(strcmp(recformat(i).type,{ALLOWED_TYPES(:).type}))
    error('recformat contains unsupported data type');
  end
  if ~isnumeric(recformat(i).count) || ~isreal(recformat(i).count) || ...
      ~isscalar(recformat(i).count) || ~isfinite(recformat(i).count) || ...
      ~(recformat(i).count >= 1) || ...
      ~(recformat(i).count == round(recformat(i).count))
    error('recformat contains an invalid count');
  end
end

if ~isnumeric(nrecs) || ~isscalar(nrecs) || ~isreal(nrecs) || ...
  ~(round(nrecs) == nrecs) || ~(nrecs >= 1)
  error('nrecs must be a positive integer');
end

finfo = dir(filename);
try
  if(finfo.bytes == 0)
    error('file %s is empty',filename);
  end
catch
  error('file %s can not be found',filename);
end
try
  [f, osmsg] = fopen(filename,'r',machineformat);
catch
  error('%s is not a valid enddianness code',machineformat);
end
if (f == -1)
  error('could not open file %s: %s',filename,osmsg);
end

% Find the end of the file header
if isempty(eoh)
  eoh_pos = 0;
else
  % If the end-of-header character sequence is not within the first
  % MAX_HEADER_SIZE bytes, there is probably a problem. Abort immediately;
  % don't scan through the remainder of the file to find the header, because
  % the user will be frustrated by having to wait so long only to receive an
  % error message.
  MAX_HEADER_SIZE = 1e4;
  % The second output returned by REGEXP is end_idx of the match
  [junk, eoh_pos] = regexp(fread(f,MAX_HEADER_SIZE,'*char')',eoh);
  if (numel(eoh_pos) > 1)
    error('more than one match was found for EOH regular expression %s',eoh);
  elseif isempty(eoh_pos)
    error('could not find EOH regular expression %s in first %d bytes', ...
        eoh,MAX_HEADER_SIZE);
  end
  frewind(f);
end

% initialize recstruct with empty fields
nfields = length(recformat);
for i = 1:nfields
  recstruct.(recformat(i).name) = zeros([0 recformat(i).count], ...
      recformat(i).type);
end

% if the file consists of nothing but a header, return the empty initialized
% recstruct (to check this, we seek to the end of the file and compare current
% position to eoh_pos)
fseek(f,0,'eof');
if (ftell(f) == eoh_pos)
  warning('file %s contains zero records',filename);
  return;
end
frewind(f);

% store the number of *complete* records that were successfully read, which may
% turn out to be less than the specified nrecs (e.g., when nrecs == Inf)
count = nrecs; 
for i = 1:nfields
  % Seek to ith field in the first record. Note that SIZE_OF is defined as an
  % inline function above.
  if (i == 1)
    start_pos = eoh_pos;
  else
    start_pos = eoh_pos + sum(cellfun(SIZE_OF,{recformat(1:(i-1)).type}));
  end
  if (fseek(f,start_pos,'bof') == -1)
    error('could not seek to position %d: %s',start_pos,ferror(f));
  end
  % Construct arguments to pass to FREAD
  fread_precision = sprintf('%d*%s=>%s', ...
      recformat(i).count,recformat(i).type,recformat(i).type);
  fread_size = [recformat(i).count nrecs];
  % To compute fread_skip, add up byte sizes of the other fields in the record.
  not_i = find(1:nfields ~= i);
  if isempty(not_i)
    fread_skip = 0;
  else
    fread_skip = sum([recformat(not_i).count] .* ...
        cellfun(SIZE_OF,{recformat(not_i).type}));
  end
  % read the ith field and add it to recstruct
  try
    [tmpdata, tmpcount] = fread(f, ...
        fread_size,fread_precision,fread_skip,machineformat);
    % Update count if we encountered an incomplete record; note that 
    % FREAD returns the total number of data elements read in tmpcount and pads
    % tmpdata with zeros for missing elements. We need to divide tmpcount by the
    % expected number of data elements in the field to get the number of records
    % whose fields were read completely
    count = min(count,floor(tmpcount/recformat(i).count));
    if (count == 0)
      error('could not read data field from record');
    end
    % We use MATLAB's dynamic field name syntax for referencing struct.
    % Transpose because FREAD returns stuff in an inconvenient orientation.
    if (size(tmpdata,2) ~= recformat(i).count)
      recstruct.(recformat(i).name) = tmpdata';
    else
      recstruct.(recformat(i).name) = tmpdata;
    end
  catch
    error('could not read data from file: %s',ferror(f));
  end
end

% If the last record in the file was incomplete, fields of recstruct may not
% match in size. Trim data fields so that they all contain COUNT rows
for i = 1:nfields
  tmp = recstruct.(recformat(i).name);
  if (size(tmp,1) > count)
    recstruct.(recformat(i).name) = tmp(1:count,:);
  end
end

% spit out header if requested
if (nargout==3)
  frewind(f);
  header = fread(f,eoh_pos,'*char')';
end
  
if (fclose(f) == -1)
  error('could not close file %s: %s',filename,ferror(f));
end

