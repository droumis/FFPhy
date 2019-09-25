% PROCOPTIONS - Processes options passed to a Matlab function.
%                   This function provides a simple means of
%                   parsing attribute-value options. Options are
%                   extracted from the namespace of the calling
%                   function.  Thus, to declare a new option, simply
%                   declare a variable in the calling function before
%                   calling the procOptions() function.  The default
%                   value of the variable (its initial value) will be
%                   changed if the arguments include a string/paramter
%                   pair in which the string corresponds to the
%                   variable name.
%
%                   This scheme was first shown to me in the code of
%                   Maneesh Sahani.  This present code was written by
%                   Caleb Kemere, 2004.
%
% Usage:  
%           function [outputs] = funfun(inputs, varargin)
%
%           var1 = default_var1;
%           var2 = default_var2;
%             ...
%           varN = default_varN;
%           unused= procOptions(varargin)
%
% Arguments:   
%            args            - a cell array of input arguments, such
%                              as that provided by VARARGIN.  Its contents
%                              should alternate between strings and
%                              values.
% Returns:
%            unused          - an optional cell array of those 
%                              string-value pairs that were unused;
%                              if this is not supplied, then a
%                              warning will be issued for each
%                              option in args that lacked a match.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [unused] = procOptions(args, varargin)

if length(args) == 0
  unused = [];
  return;
end

warn = 0;
if (length(varargin) >= 2)
   if strcmpi(varargin{1},'warn')
      warn = str2num(varargin{2});
   end
end
vars = evalin('caller','whos');

% Check the number of input arguments
n = length(args);
if (mod(n, 2))
  error('Each option must be a string/value pair.');
end

% Now process all arguments
nunused = 0;
unused = {};
i = 1;
while (i <= length(args))
  found = 0;
  for j=1:length(vars)
    if strcmpi(args{i}, vars(j).name)
      assignin('caller', vars(j).name, args{i + 1});
      found = 1;
      break;
    end
  end
  if (~found)
    if (warn)
      warning(sprintf('Option ''%s'' not used.', args{i}));
      args{i}
    end
    nunused = nunused + 1;
    unused{2*nunused-1} = args{i};
    unused{2 * nunused} = args{i + 1};
  end
  i = i + 2;
end

