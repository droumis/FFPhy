function varargout = subsref(t,s)
%SUBSREF Subscripted reference for a piecewisedistribution object.
%   B = SUBSREF(T,S) is called for the syntax T(X) when T is a
%   piecewisedistribution object.  S is a structure array with the fields:
%       type -- string containing '()', '{}', or '.' specifying the
%               subscript type.
%       subs -- Cell array or string containing the actual subscripts.
%
%   P = T(X), where T is a piecewisedistribution object, returns the cdf
%   of the piecewise distribution evaluated at the values in X.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2006/12/15 19:32:54 $

switch s(1).type
case '()'
    if ~isscalar(s)
        error('stats:stattree:subsref:InvalidSubscriptExpr', ...
              'Subscript expression is invalid.  Use () indexing to evaluate the cdf.');
    
    elseif numel(s(1).subs) ~= 1
        error('stats:stattree:subsref:BadSubscript', ...
              '() indexing requires a single array of input values.');
    end
       
    % '()' is a reference to the cdf method
    [varargout{1:nargout}] = cdf(t,s(1).subs{1});


case '{}'
    error('stats:piecewisedistribution:subsref:CellSubscript', ...
          'The PIECEWISEDISTRIBUTION class does not support cell array indexing.');
        
case '.'
    % Support dot subscripting to call methods
    if isscalar(s)
        args = {};
    elseif numel(s)>2 || ~isequal(s(2).type,'()')
         error('stats:piecewisedistribution:subsref:BadSubscript', ...
               'Invalid subscripting for PIECEWISEDISTRIBUTION object.');
    else   % numel(s)==2
        args = s(2).subs;
    end
    
    % For improved performances these methods are hard-coded instead of
    % being found via the "methods" function and called via "feval"
    switch(s(1).subs)
      case 'cdf'
        [varargout{1:nargout}] = cdf(t,args{:});
      case 'icdf'
        [varargout{1:nargout}] = icdf(t,args{:});
      case 'random'
        [varargout{1:nargout}] = random(t,args{:});
      case 'pdf'
        [varargout{1:nargout}] = pdf(t,args{:});
      case 'segment'
        [varargout{1:nargout}] = segment(t,args{:});
      case 'nsegments'
        [varargout{1:nargout}] = nsegments(t,args{:});
      case 'boundary'
        [varargout{1:nargout}] = boundary(t,args{:});
      case 'disp'
        [varargout{1:nargout}] = disp(t,args{:});
        
      otherwise
         error('stats:piecewisedistribution:subsref:BadMethod', ...
               'Invalid method name "%s".',s(1).subs);
    end
end
