% function h = seplot(x, varargin)
% 
% Copyright 2004 Caleb Kemere
% May be freely distributed.
%
% Either called with
% seplot(x, y)
% or
% seplot(y)
% 
% Plots the mean of y vs. x with confidence bounds.
% 
% Optional Arguments:
%
% 'lineStyle' - use a matlab linestyle such as '--' or ':' (default: '-')
% 'lineWidth' - width of the main line (default: 1)
% 'marker'    - markers for the points on the main line (default: '')
% 'color'     - color of the main line (default: matlab default, i.e., 'b') 
%                (note that the confidence bounds are fcolor)
%
% 'varType'   - the type of confidence bound - could be:
%                 'se'   - standard error of the mean
%                 'std'  - standard deviation
%                 'conf' - percentiles (spec'd by 'confAlpha')
% 'confAlpha' - the percentiles to use if varType is 'conf' (default: 0.05)
% 'fcolor'    - color to use for confidence intervals (default: faded value of 'color')
% 'faceAlpha' - alpha for confidence intervals (default: [])
%                (note that alpha blending doesn't export nicely
%                 in vector formats, so this is good for display
%                 or png/jpeg exports, but this must be recreated
%                 post-export for vector files)
%
% Excess arguments are passed on to the 'plot' command for the mean.
% 
%

function h = seplot(x, varargin)
if ( (nargin > 1) & ~ischar(varargin{1}) )
   y = varargin{1};
   hh = plot(x,mean(y));
   flag2 = 1;
   if (length(varargin) > 1)
      varargin = {varargin{2:end}};
   else
      varargin = {};
   end
else
   y = x;
   hh = plot(mean(y));
   flag2 = 0;
end


holdState = ishold;
hold on
[lineStyle, lineWidth, color, marker, confAlpha, faceAlpha, fcolor, vartype, otherArgs] = ...
process_options(varargin, ...
   'lineStyle', '-', 'lineWidth', 1, 'Color', [], 'marker', [], ...
   'confAlpha', 0.05, 'faceAlpha', [], 'fcolor', [], 'VarType', 'se');


if isempty(vartype) | strcmpi(vartype,'se')
   se = transpose(std(y)./sqrt(size(y,1) - 1));
elseif strcmpi(vartype,'std')
   se = transpose(std(y));
elseif strcmpi(vartype,'conf')
   se = [];
   se1 = prctile(y,100*(1-confAlpha/2));
   se2 = prctile(y,100*confAlpha/2);
end

xdata = get(hh,'XData');
ydata = get(hh,'YData')';
longXdata = [xdata, fliplr(xdata)];
if ~isempty(se)
 longSE = [ydata + se; flipud(ydata - se)];
else
 longSE = [se1(:); flipud(se2(:))];
end
if isempty(fcolor)
 fcolor = get(hh,'Color');
end
ff = fill(longXdata,longSE,fcolor,'EdgeColor','none');
delete(hh);

my = mean(y);

if flag2
   if (length(otherArgs) > 0)
      hh = plot(x,my,'linestyle',lineStyle,'linewidth',lineWidth, ...
         otherArgs{:});
   else
      hh = plot(x,my,'linestyle',lineStyle,'linewidth',lineWidth);
   end
else
   if (length(otherArgs) > 0)
      hh = plot(my,'linestyle',lineStyle,'linewidth',lineWidth, ...
         otherArgs{:});
   else
      hh = plot(my,'linestyle',lineStyle,'linewidth',lineWidth);
   end
end

if ~isempty(color)
  set(hh,'Color',color);
  hsl_color = rgb2hsl(color);
  hsl_color(3) = min(1,hsl_color(3) + 0.5);
  faded_color = hsl2rgb(hsl_color);
  set(ff,'FaceColor',faded_color);
else
  set(hh,'Color',fcolor*0.5);
end
if ~isempty(marker)
  set(hh,'Marker',marker);
end
if ~isempty(faceAlpha)
  set(ff,'FaceAlpha',faceAlpha);
end

h = [hh; ff'];

if (holdState == 0)
   hold off
end


% PROCESS_OPTIONS - Processes options passed to a Matlab function.
%                   This function provides a simple means of
%                   parsing attribute-value options.  Each option is
%                   named by a unique string and is given a default
%                   value.
%
% Usage:  [var1, var2, ..., varn[, unused]] = ...
%           process_options(args, ...
%                           str1, def1, str2, def2, ..., strn, defn)
%
% Arguments:   
%            args            - a cell array of input arguments, such
%                              as that provided by VARARGIN.  Its contents
%                              should alternate between strings and
%                              values.
%            str1, ..., strn - Strings that are associated with a 
%                              particular variable
%            def1, ..., defn - Default values returned if no option
%                              is supplied
%
% Returns:
%            var1, ..., varn - values to be assigned to variables
%            unused          - an optional cell array of those 
%                              string-value pairs that were unused;
%                              if this is not supplied, then a
%                              warning will be issued for each
%                              option in args that lacked a match.
%
% Examples:
%
% Suppose we wish to define a Matlab function 'func' that has
% required parameters x and y, and optional arguments 'u' and 'v'.
% With the definition
%
%   function y = func(x, y, varargin)
%
%     [u, v] = process_options(varargin, 'u', 0, 'v', 1);
%
% calling func(0, 1, 'v', 2) will assign 0 to x, 1 to y, 0 to u, and 2
% to v.  The parameter names are insensitive to case; calling 
% func(0, 1, 'V', 2) has the same effect.  The function call
% 
%   func(0, 1, 'u', 5, 'z', 2);
%
% will result in u having the value 5 and v having value 1, but
% will issue a warning that the 'z' option has not been used.  On
% the other hand, if func is defined as
%
%   function y = func(x, y, varargin)
%
%     [u, v, unused_args] = process_options(varargin, 'u', 0, 'v', 1);
%
% then the call func(0, 1, 'u', 5, 'z', 2) will yield no warning,
% and unused_args will have the value {'z', 2}.  This behaviour is
% useful for functions with options that invoke other functions
% with options; all options can be passed to the outer function and
% its unprocessed arguments can be passed to the inner function.

% Copyright (C) 2002 Mark A. Paskin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = process_options(args, varargin)

% Check the number of input arguments
n = length(varargin);
if (mod(n, 2))
  error('Each option must be a string/value pair.');
end

% Check the number of supplied output arguments
if (nargout < (n / 2))
  error('Insufficient number of output arguments given');
elseif (nargout == (n / 2))
  warn = 1;
  nout = n / 2;
else
  warn = 0;
  nout = n / 2 + 1;
end

% Set outputs to be defaults
varargout = cell(1, nout);
for i=2:2:n
  varargout{i/2} = varargin{i};
end

% Now process all arguments
nunused = 0;
i = 1;
while (i < length(args))
  found = 0;
  for j=1:2:n
    if strcmpi(args{i}, varargin{j})
      varargout{(j + 1)/2} = args{i + 1};
      found = 1;
      break;
    end
  end
  if (~found)
    if (warn)
      warning(sprintf('Option ''%s'' not used.', args{i}));
    else
      nunused = nunused + 1;
      unused{nunused} = args{i};
      % unused{2*nunused-1} = args{i};
      % unused{2 * nunused} = args{i + 1};
    end
    if (i+1 == length(args))
      nunused = nunused + 1;
      unused{nunused} = args{i+1};
    end
    i = i - 1;
  end
  i = i + 2;
end

% Assign the unused arguments
if (~warn)
  if (nunused)
    varargout{nout} = unused;
  else
    varargout{nout} = cell(0);
  end
end
