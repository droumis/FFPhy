function nice
%NICE              Custom plot prettification.
%   The function NICE is a macro for Matlab graphics settings, and its
%   behavior depends on the contents of the current axes.
%
%   If the plot contains only 2-D (i.e., empty ZData) line objects NICE
%   sets the Y-limits depending on the data:
%               Data Range        Gives Limits
%                [ 0,1]           [-0.2 1.2]
%                [-1,1]           [-1.2 1.2]
%                [ 0,B]           [-0.1 1.2]*B
%                [-A,B]           [-1.2 1.2]*C, where C=MAX(ABS(A,B))
%   where the first condition that matches is used.
%
%   If the plot contains surface/patch objects, 3-D visualization aids are
%   set, including: two lights, interpolated shading (if possible), and
%   constant view angle for 3-D rotation.


%%%%%%%%%%%%%% Actually, it doesn't do the following for now ... not sure what the best thing is
%       The X-limits are set similarly if the X-values are not uniformly spaced;
%       Uniformly spaced x-values cause the X-limits to be set between the min
%       X-value - 1 and the max X-value + 1.

childs = get(gca, 'Children');
if (length(childs) == 1),  childs = childs([1,1]);  end;        % ensure get(childs, prop) always returns cell
childs(ismember(get(childs, 'Type'), {'text', 'light'})) = [];  % ignore these
childtypes = get(childs, 'Type');

if (all(strcmp(childtypes, 'line')))
	
	twoD = cellfun('isempty', get(childs, {'ZData'}));
	if (all(twoD))
		ydata = get(childs, 'YData');   ydata = cat(2, ydata{:});
		xdata = get(childs, 'XData');   xdata = cat(2, xdata{:});
		
		miny = min(ydata);  maxy = max(ydata);
		maxabs    = max(abs([miny,maxy,1]));
		if (miny < 0),   ylim = [-1.2, 1.2]*maxabs;
		else,            ylim = [-0.2, 1.2]*maxabs;
		end

		xlim = [min(xdata)-1, max(xdata)+1];   % temporary -- just tightens the X limits.
		
		set(gca, 'YLim', ylim, 'XLim', xlim);
		zoom reset;
	end
	
elseif (any(ismember(childtypes, {'surface', 'patch'})))
	childs(strcmp(childtypes, 'line')) = [];
	
	axis tight;  set(gca,'CameraViewAngleMode','manual', 'CameraViewAngle', 9);

	interpable = ~cellfun('isempty', get(childs, {'CData'}));
	set(childs(interpable), 'FaceColor', 'interp', 'EdgeColor', 'none');

	material dull;	
	lighting gouraud;   % not as nice as phong lighting, but faster to render
	lightangle(90,60);  lightangle(270,10);
	
end