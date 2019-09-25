function [ ix ] = peakpick(x,choice);
% PEAKPICK peak picking.
% PEAKPICK(X) picks the peaks of a sequence.

if nargin < 2
   choice = 'max'; % default choice
end
switch lower(choice)
   case {'max'}
      ix = find(sign(-sign(diff(sign(diff(x)))+0.5)+1))+1;
   case {'min'}
      ix = find(sign(-sign(diff(sign(diff(-x)))+0.5)+1))+1;
   otherwise
      warning('Default search of local maximum !');
      ix = find(sign(-sign(diff(sign(diff(x)))+0.5)+1))+1;
end;


end

