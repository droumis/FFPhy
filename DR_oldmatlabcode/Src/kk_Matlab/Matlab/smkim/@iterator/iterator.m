classdef iterator < handle
%ITERATOR Class for accessing struct array data that is distributed across multiple files
%
%   IT = ITERATOR(FILENAME_PATTERN,VARIABLE_NAME,VALIDATION_FUNC) returns an
%   iterator object which can be used to access column-vector struct array data.
%   FILENAME_PATTERN must be a string or a cell array of strings which can be
%   expanded into filenames by the DIR function. Each of the target files must
%   be a .MAT file which, when loaded into the workspace, produces a
%   column-vector (Nx1) struct array variable named VARIABLE_NAME. For type
%   validation, VALIDATION_FUNC(VARIABLE_NAME) must return true. Furthermore,
%   all of the struct arrays that are loaded from these files must have the
%   exact same fields (no missing or extra fields).
%
%   If the iterator constructor is successful, then the data can then be
%   accessed as if they belonged to single concatenated column struct array.
%   Examples:
%     {IT.fieldname} returns a cell array of the contents of fieldname for all
%     of the struct elements in the dataset.
%     IT(end-2:end) returns the last three struct elements in the dataset.
%     IT([IT.fieldname] > 2) returns all elements for which the value in
%     fieldname is greater than 2 (assuming that the GT comparison function
%     returns successfully)
%
%   Instead of ARRAYFUN, use the MAP method:
%     MAP(IT,@my_function)
%
%   The FIND method is overloaded:
%     FIND(IT,@test_condition) returns the indices of IT for which
%     @test_condition evaluates to true.
%
%   Performance may be slow if many files have to be loaded. Also, vectorized
%   indexing access may fail if a memory limit is exceeded.
%
%   Known anomalous behavior: IT(index_vector).fieldname does not return a
%   comma-separated list, but instead returns a cell array. This works around a
%   bug in the way that MATLAB determines nargout when the subsref method is
%   overloaded. 
%
%Written by SMK, 2009 November 10.
%

  properties (Access = private, Constant)
    % How many bytes can safely be returned when calling subsref? This constant
    % should be modified according to the platform.
    max_argout_bytes = 1e9; 
    % Which fields do we want to store in self.lookup for fast access?
    allowed_metadata_fields = { ...
        'subject'                 ; ...
        'day'                     ; ...
        'epoch'                   ; ...
        'environment'             ; ...
        'timerange'               ; ...
        'depth'                   ; ...
        'hemisphere'              ; ...
        'region'                  ; ...
        'Fs'                      ; ... 
        'reference'               ; ... 
        'electrode'               ; ...
        'tetrode'                 ; ...
        'uid'                     };
  end

  properties (SetAccess = private, GetAccess = public)
    % pattern used to match filenames for loading
    filename_pattern
    % name of the workspace variable that is created when the files are loaded
    variable_name
    % function used to validate the loaded workspace variable
    validation_func
    % fieldnames of the underlying struct arrays
    fields
    % fieldnames that are included in metadata
    metadata_fields
    % filenames is a cell array of the unique .filename fields of lookup
    filenames
  end

  properties (Access = private)
    % struct for holding temporary copy of data struct
    tmpdata
    % filename that was loaded as tmpdata
    last_loaded
    % lookup(i).filename is the filename to be loaded for the ith item of the
    % iterator object, and lookup(i).index is the index of the struct element
    % in the struct array that is loaded from filename which corresponds to the
    % ith item
    lookup
    % a struct array that contains only those fields that are specified in
    % self.allowed_metadata_fields, allowing fast access to these fields without
    % having to load the entire *.mat files
    metadata
    % indices is a cell array of vectors of indices, such that indices{i} lists
    % the elements of lookup(i) whose filename field matches filenames{i}
    indices
  end

  methods (Access = public)

    function self = iterator(filename_pattern,variable_name,validation_func)
      if (nargin ~= 3)
        error('iterator constructor must be called with three arguments');
      end
      if ischar(filename_pattern)
        filename_pattern = {filename_pattern};
      end
      if ~iscellstr(filename_pattern) && ...
          ~all(cellfun(@(c) isdir(fileparts(c)),filename_pattern))
        error('filename_pattern must contain fully-qualified filesystem path');
      end
      self.filename_pattern = filename_pattern;
      if ~ischar(variable_name)
        error('variable_name must be a string');
      end
      self.variable_name = variable_name;
      if ~isa(validation_func,'function_handle') || ...
          (nargin(validation_func) ~= 1)
        error(['validation_func must be a function handle that takes a ' ...
            'single input'])
      end
      self.validation_func = validation_func;
      self.fields = [];
      self.metadata_fields = [];
      self.lookup = struct('filename',{},'index',{});
      self.metadata = [];
      for i = 1:numel(filename_pattern)
        files_list = dir(filename_pattern{i});
        if (numel(files_list) < 1)
          warning('no files that match pattern %s were found', ...
              filename_pattern{i});
        end
        for j = 1:numel(files_list)
          filename = sprintf('%s/%s', ...
              fileparts(filename_pattern{i}),files_list(j).name);
          try
            % this has side effects of assigning self.tmpdata and
            % self.last_loaded
            self.load_data(filename);
          catch
            error(['could not load expected data struct named %s from ' ...
                'candidate file %s'],variable_name,filename);
          end
          if isequal(self.fields,[])
            self.fields = fieldnames(self.tmpdata);
            % select fields for inclusion in self.metadata
            self.metadata_fields = intersect(self.allowed_metadata_fields, ...
                self.fields);
            tmpargs = repmat({{}},[1 2*numel(self.metadata_fields)]);
            [tmpargs{1:2:end}] = deal(self.metadata_fields{:});
            self.metadata = struct(tmpargs{:});
          end
          for k = 1:numel(self.tmpdata)
            self.lookup(end+1,1) = struct('filename',filename,'index',k);
            self.metadata(end+1,1) = orderfields(rmfield( ...
                self.tmpdata(k),setdiff(self.fields,self.metadata_fields)), ...
                self.metadata_fields);
          end
        end
      end
      assert(isequal(size(self.lookup),size(self.metadata)));
      if isempty(self.lookup)
        error('iterator is empty');
      end
      % for each unique filename, find matching entries in self.lookup
      self.filenames = unique({self.lookup.filename});
      self.indices = cellfun(@(c) find(strcmp(c,{self.lookup.filename})), ...
          self.filenames,'UniformOutput',false);
      assert(isequal(size(self.filenames),size(self.indices)));
    end

    % The subsref method is defined in a separate m-file, subsref.m
    varargout = subsref(self,s) 

    % The subsasgn method is disabled
    function varargout = subsasgn(self,varargin)
      error('iterator class is read-only, so the subsasgn method is disabled');
    end

    % Correct behavior of numel method is very important! MATLAB documentation:
    % "In the case of the overloaded subsref function for brace and dot
    % indexing ... numel is used to compute the number of expected outputs
    % (nargout) returned from subsref."
    function val = numel(self,varargin)
      if (nargin == 1)
        val = numel(self.lookup);
      elseif (nargin > 1)
        val = prod(cellfun(@numel,varargin));
      end
    end

    % end method depends on the numel method
    function val = end(self,k,n)
      if (k == 1)
        val = numel(self);
      else
        val = 1;
      end  
    end

    % size method depends on the numel method
    function varargout = size(self,varargin)
      if (nargin == 1)
        if nargout <= 1
          varargout{1} = [numel(self), 1];
        else
          varargout{1} = numel(self);
          [varargout{2:nargout}] = deal(1);
        end
      elseif (nargin == 2)
        dim = varargin{1};
        if (dim == 1)
          varargout{1} = numel(self);
        else
          % all dimensions beyond the first are singleton dimensions
          varargout{1} = 1;
        end
      else
        error('Too many input arguments');
      end
    end

    % length method depends on the size method
    function val = length(self)
      val = max(size(self));
    end

    % ndims method depends on the size method
    function val = ndims(self)
      val = length(size(self));
    end
   
    % fieldnames method requires self.fields to be correct
    function val = fieldnames(self)
      val = self.fields;
    end

    function val = isfield(self, field)
      if ischar(field)
        if any(strcmp(field,self.fields))
          val = true;
        else
          val = false;
        end
      elseif iscellstr(field)
        val = false(size(field));
        for i = 1:numel(field)
          if any(strcmp(field{i},self.fields))
            val(i) = true;
          end
        end
      else
        % This may seem weird but is consistent with the behavior of the builtin
        % isfield function
        val = false;
      end
    end

    % This method returns a struct array consisting of elements that satisfy a
    % test condition
    function val = find(self, test_condition)
      if ~isa(test_condition,'function_handle');
        error('test_condition must be a function handle');
      end
      try
        % call map method
        tf = self.map(test_condition);
        assert(islogical(tf) && isvector(tf) && isequal(size(tf),size(self)));
      catch
        error(['test_condition failed to return logical vector of the same ' ...
            'size as %s'],inputname(1));
      end
      val = find(tf);
    end

    % This method acts as a wrapper around arrayfun, taking the same optional
    % parameter/value pairs
    function val = map(self,func,varargin)
      if ~isa(func,'function_handle');
        error('func must be a function handle');
      end
      p = inputParser();
      p.addParamValue('UniformOutput',true,@(v) islogical(v) && isscalar(v));
      p.addParamValue('ErrorHandler', ...
          @(varargin) error('Error in map method'), ...
          @(v) isa(v,'function_handle') && ...
          ((nargin(v) == 2) || (nargin(v) < 0)));
      try
        p.parse(varargin{:});
      catch
        error('Invalid parameter/value pair(s)');
      end
      % Attempt to evaluate the function on metadata first, and if this fails
      % then iterate over files
      try
        val = arrayfun(func,self.metadata, ...
            'UniformOutput',p.Results.UniformOutput);
      catch
        for i = 1:numel(self.filenames)
          self.load_data(self.filenames{i});
          try
            if p.Results.UniformOutput
              val(self.indices{i},1) = arrayfun(func, ...
                  self.tmpdata([self.lookup(self.indices{i}).index]), ...
                  varargin{:});
            else
              out = arrayfun(func, ...
                  self.tmpdata([self.lookup(self.indices{i}).index]), ...
                  varargin{:});
              [val{self.indices{i}}] = out{:};
            end
          catch exception
            rethrow(exception);
          end
        end
      end
    end

    % This method returns elements that satisfy a test condition
    function val = select(self,test_condition)
      if ~isa(test_condition,'function_handle');
        error('test_condition must be a function handle');
      end
      try
        tf = map(self,test_condition);
        assert(islogical(tf) && isvector(tf) && (numel(tf) == numel(self)));
      catch
        error('test_condition does not evaluate to a logical vector');
      end
      val = subsref(self,struct('type',{'()'},'subs',{{tf}}));
    end

    function disp(self)
      disp(sprintf('%dx1 %s iterator with fields:', ...
          numel(self),self.variable_name));
      disp(sprintf('    %s\n',self.fields{:}));
    end

    function display(self)
      self.disp();
    end

    % The get method allows access to properties without relying on dot
    % notation, which can be ambiguous when a data field has the same name as a
    % property
    function val = get(self,propname)
      if ismember(propname,properties(self))
        val = self.(propname);
      end
    end

  end

  methods (Access = private)

    % Wrapper method for loading data from file, maintains state of
    % self.last_loaded and self.tmpdata
    function load_data(self,filename)    
      if ~ischar(filename) || ~isdir(fileparts(filename)) || ...
          (exist(filename) ~= 2)
        error('%s is not a valid filename');
      end
      if strcmp(self.last_loaded,filename)
        % if filename matches self.last_loaded, do nothing; data is already
        % loaded as self.tmpdata
        return;
      end
      load(filename);
      if (exist(self.variable_name) ~= 1)
        error('file %s does not load a variable called %s', ...
            filename,self.variable_name);
      end
      if ~isstruct(eval(self.variable_name)) || ...
          ~isvector(eval(self.variable_name)) || ...
          (size(eval(self.variable_name),2) > 1)
        error('%s loaded from %s must be a column-vector struct array', ...
            self.variable_name,filename);
      end
      if ~self.validation_func(eval(self.variable_name))
        error('%s loaded from %s is not a valid %s data struct', ...
            self.variable_name,filename,self.variable_name);
      end
      % Reorder fields to allow comparisons
      if ~isequal(self.fields,[])
        try
          self.tmpdata = orderfields(eval(self.variable_name),self.fields);
          assert(isequal(fieldnames(self.tmpdata),self.fields))
        catch
          error(['%s data struct loaded from %s is missing fields ' ...
              'or has unexpected extra fields'],self.variable_name,filename);
        end
      else
        self.tmpdata = eval(self.variable_name);
      end
      self.last_loaded = filename;
    end

  end

end



