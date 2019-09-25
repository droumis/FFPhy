function h = fig(id, varargin)
% fig - find or create a figure
%
% 	FIG and FIG(H) where H is a figure handle are identical to 
%	the built-in FIGURE.  
%
%	FIG('name') looks for a figure with the specified name.  If it
%	finds one it makes it current; if it fails it creates one with
%	the given name.
%
%	FIG(ID, 'property', 'value', ...) obtains a handle to the
%	figure with the specified ID (either handle or name), creating
%	it if necessary, and sets the specfied properties.

if nargin < 1
  h = figure;
  return;
end

if (~ischar(id))
  figure(id);
  h = id;
else
  h = findobj ('type', 'figure', 'name', id);
  if (~ isempty(h))
    figure(h);
  else
    h = figure('name', id);
  end
end

set(h, varargin{:});
