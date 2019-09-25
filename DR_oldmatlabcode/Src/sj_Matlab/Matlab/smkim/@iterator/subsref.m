function varargout = subsref(self,s)
%SUBSREF Overloaded subsref method for the iterator class
%
%Known bug: 
% This implementation is complicated by the fact that MATLAB does not assign the
% correct nargout value when subsref is implicitly called by
% iterator_obj(indices).fieldname; nargout should equal numel(indices) but
% instead it is zero. As a workaround, in these cases the return value (which
% should be a comma-separated list) is packaged in a cell array.
% 

  switch s(1).type

  case '{}'
    error('Cell contents reference from a non-cell array object');

  case '.'
    % If s(1).type is the dot reference '.', then s(1).subs must be a
    % an element of self.fields. The .method() call syntax and .property get
    % syntax are suppressed to avoid namespace collisions with the data fields.
    field = s(1).subs;
    if ismember(field,properties(self))
      warning(['Access to iterator property ''%s'' by dot notation is ' ...
          'disabled to avoid namespace collision with data field names. ' ...
          'Try using the get method get(%s,''%s'') instead.'], ...
          field,inputname(1),field);
    end
    if ismember(field,methods(self))
      warning(['Calling iterator method ''%s'' by dot notation is ' ...
          'disabled to avoid namespace collision with data field names. ' ...
          'Try calling the method as %s(%s,...) instead.'], ...
          field,field,inputname(1));
    end
    if ~ismember(field,self.fields)
      error('Reference to non-existent field ''%s''',field);
    end
    % Field reference is processed differently if iterator is a scalar
    if (numel(self) == 1)
      self.load_data(self.lookup.filename);
      try
        varargout{1} = subsref(self.tmpdata(self.lookup.index),s);
      catch exception
        rethrow(exception);
      end
    else
      if (length(s) > 1)
        error(['Field reference for multiple structure elements that ' ...
            'is followed by more reference blocks is an error']);
      end
      if ismember(field,self.metadata_fields)
        % We can lookup using self.metadata, avoiding costly load() calls
        [varargout{1:numel(self)}] = self.metadata.(field);
      else
        for i = 1:numel(self.filenames)
          self.load_data(self.filenames{i});
          [varargout{self.indices{i}}] = ...
              self.tmpdata([self.lookup(self.indices{i}).index]).(field);
          meminfo = whos('varargout');
          assert((numel(meminfo) == 1) && ...
              strcmp('iterator.subsref',meminfo(1).nesting.function));
          if (meminfo(1).bytes > self.max_argout_bytes)
            error(['Exceeded memory limit (%d bytes). Try iterating over ' ...
                'elements one-by-one inside of a for loop instead of ' ...
                'vector indexing'],self.max_argout_bytes);
          end
        end
      end
    end
    
  case '()'
    % If s(1).type is an array index reference '()', then s(1).subs must be
    % convertible to linear indices in the range (1:numel(self)), and (length(s)
    % == 1) || strcmp('.',s(2).type)
    if (length(s) > 1) && ~strcmp('.',s(2).type)
      error('Error: ()-indexing must appear last in an index expression');
    end    
    ind = s(1).subs;
    for dim = 1:length(ind)
      % Convert logical indices
      if islogical(ind{dim})
        ind{dim} = find(ind{dim});
      end
      % Expand any colon subscripts
      if ischar(ind{dim}) && strcmp(':',ind{dim})
        ind{dim} = 1:size(self,dim);
      end
    end
    if any(cellfun(@(i) any(i <= 0) || ~isreal(i) || ~all(isfinite(i)) || ...
        ~all(round(i) == i),ind))
      error(['Subscript indices must be either real positive integers or ' ...
          'logicals']);
    end
    if any(arrayfun(@(dim) ~isempty(ind{dim}) && ...
        (max(ind{dim}) > size(self,dim)),1:length(ind)))
      error('Index exceeds object dimensions; size(%s)=%s', ...
          inputname(1),mat2str(size(self)));
    end
    % Convert to linear index. Because the higher dimensions are all singletons,
    % we are guaranteed that the index along the first dimension is equivalent
    % to a linear index. Elements of ind map one-to-one to elements of
    % varargout.
    ind = ind{1};
    switch (numel(ind))
    case 0
      if (length(s) == 1)
        % Return a 0x0 struct array with the correct fields
        tmpargs = repmat({{}},[1 2*numel(self.fields)]);
        [tmpargs{1:2:end}] = deal(self.fields{:});
        varargout{1} = struct(tmpargs{:});
      else
        % This is consistent with MATLAB builtin behavior
        varargout = {};
      end
    case 1
      % If ind is a single element, then we grab a scalar struct and pass along
      % any remaining subscript references
      self.load_data(self.lookup(ind).filename);
      if (length(s) > 1)
        try
          varargout{1} = subsref(self.tmpdata(self.lookup(ind).index),s(2:end));
        catch exception
          rethrow(exception);
        end
      else
        varargout{1} = self.tmpdata(self.lookup(ind).index);
      end
    otherwise
      if (length(s) > 1) && (any(strcmp('()',{s(2:end).type})) || ...
          any(strcmp('{}',{s(2:end).type})))
        error('Error: ()-indexing must appear last in an index expression');
      end
      switch (length(s))
      case 2
        field = s(2).subs;
        if ~ismember(field,self.fields)
          error('Reference to non-existent field ''%s''',field);
        end
        if ismember(field,self.metadata_fields)
          % Satisfy this query using self.metadata
          [varargout{1:numel(ind)}] = self.metadata(ind).(field);
        else
          filenames = unique({self.lookup(ind).filename});
          for i = 1:numel(filenames)
            self.load_data(filenames{i});
            out_idx = find(strcmp(filenames{i},{self.lookup(ind).filename}));
            try
              [varargout{out_idx}] = ...
                  self.tmpdata([self.lookup(ind(out_idx)).index]).(field);
            catch exception
              rethrow(exception);
            end
            meminfo = whos('varargout');
            assert((numel(meminfo) == 1) && ...
                strcmp('iterator.subsref',meminfo(1).nesting.function));
            if (meminfo(1).bytes > self.max_argout_bytes)
              error(['Exceeded memory limit (%d bytes). Try iterating over ' ...
                  'elements one-by-one inside of a for loop instead of ' ...
                  'vector indexing'],self.max_argout_bytes);
            end
          end
        end
      case 1
        % If given only ()-indexing, return a struct array of length numel(ind)
        filenames = unique({self.lookup(ind).filename});
        for i = 1:numel(filenames)
          % For each file that is loaded, populate the corresponding elements of
          % the output struct array
          self.load_data(filenames{i});
          out_idx = find(strcmp(filenames{i},{self.lookup(ind).filename}));
          outstruct(out_idx,1) = ...
              self.tmpdata([self.lookup(ind(out_idx)).index]); 
          meminfo = whos('outstruct');
          assert((numel(meminfo) == 1) && ...
              strcmp('iterator.subsref',meminfo(1).nesting.function));
          if (meminfo(1).bytes > self.max_argout_bytes)
            error(['Exceeded memory limit (%d bytes). Try iterating over ' ...
                'elements one-by-one inside of a for loop instead of ' ...
                'vector indexing'],self.max_argout_bytes);
          end
        end
        assert(numel(outstruct) == numel(ind));
        varargout{1} = outstruct;
      otherwise
        error('Scalar index required for this type of multi-level indexing'); 
      end
    end
    if (length(varargout) > 1)
      warning(['Due to a bug in the way that MATLAB determines nargout, ' ...
          'comma-separated return values are packaged in a cell array']);
      varargout{1} = varargout;
    end

  end

end % end subsref method

