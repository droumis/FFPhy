function hist3d(data_mtx, bins, smoothradius, smoothsd, levels, color)

% HIST2D  3-Dimensional Density Histogram.
% 
% HIST3D(DATA_MTX), where DATA_MTX is an (m x 3) matrix, helps to
%    visualize the density of a 3-D scatter plot.
%    The visualization is computed as follows.  The rows generate a 3-D
%    histogram which is then smoothed with a 3-D gaussian kernel.  This
%    histogram is treated as a scalar-valued volume and isocontours are
%    found at a number of density levels, each drawn as a transparent,
%    colored surface.  Higher density contours can be shaded with colors
%    higher in the current figure colormap, or all contours can be drawn
%    with a single color if so desired.  Finally, the axis tick marks
%    are colored to assist in tracking orientation (seems to work better
%    than axis labels) -- the coloring is (r,g,b) => (x,y,z).
%   *Note: REDUCEPATCHES is used before drawing the density to keep the
%    number of patches to a manageable level; this means that fine details
%    of the volume may not completely accurate.  Its just a visualization toy. 
%
%    The surface is controlled by 5 optional arguments:
%
%       HIST3D(DATA_MTX, BINS, SMOOTHRADIUS, SMOOTHSD, LEVELS, PROGRESS)
%
%    Any optional argument specified as [] will use its default value.
%    The definitions of the control arguments and there defaults are:
%      BINS            (default: 50) number of bins used to
%                       discretize the data for histogramming
%      SMOOTHRADIUS    (default: 3) linear dimension of gaussian
%                       kernel used to smooth.  Units of bins.
%      SMOOTHSD        (default: 1) standard devation of smoothing
%                       kernel.  Units of bins.
%      LEVELS          (default: 10) A scalar value is taken as the #
%                       of evenly spaced values (starting at 1 and ending
%                       below max density) at which to draw a contour.
%                       Vectors with values between [0..1] are treated
%                       as fractions of max density, e.g., [0.1 0.2 0.3]
%                       shows contours at 10%, 20% and 30% of max density.
%      COLOR           (default: evenly spaced intervals in the colormap
%                       are used for increasing density contours)  If
%                       specified as [R G B], it is taken to represent
%                       a single color at which all contours will be drawn.
%

% Process inputs and set up default values where needed.
if (size(data_mtx,2) ~= 3)
    error('Data must by [m x 3].');
end
if (nargin < 2 | isempty(bins))
    bins = 50;
end
if (nargin < 3 | isempty(smoothradius))
    smoothradius = 3;
end
if (nargin < 3 | isempty(smoothsd))
    smoothsd = 1;
end
if (nargin < 3 | isempty(levels))
    levels = 10;
else
    if (numel(levels) == 1)
        if (levels < 1)
            error('If the ''levels argument is a scalar, it must be > 1.''');
        end
    elseif ((ndims(levels) > 2) | (any(levels > 1)))
        error('If the ''levels'' argument is a vector, it must be 2-D and it cannot have any values > 1.');
    end
end

% 'patch_reduction' (not currently a control option) refers to the degree
% to which contour surfaces are to be simplified.  It defines the fraction of
% faces to be retained during the 'reducepatch' call.
patch_reduction = 0.3;

% Scale the data to a range convenient for histogramming . . .
[data1,min1,max1] = rescale_data(data_mtx(:,1),1,bins);
[data2,min2,max2] = rescale_data(data_mtx(:,2),1,bins);
[data3,min3,max3] = rescale_data(data_mtx(:,3),1,bins);

% Construct the 3-D scatter density as a smoothed 3-D histogram.
% 'sparse' is an effective histogramming tool, but it only works in 2-D so we
%  combine the first 2 dimensions to allow it to count . . .
counts = full(sparse(round(data2)+(bins)*(round(data1)-1), round(data3), 1, bins^2, bins));
counts = reshape(counts,[bins,bins,bins]);  % . . . and now separate them.
smoothdata = smooth3(counts, 'gaussian', smoothradius, smoothsd);

% Now determine which contour levels we'll want to draw.
maxdens = max(smoothdata(:));
if (numel(levels) == 1)
    showlevels = linspace(1,maxdens,levels+1);
    showlevels = showlevels(1:end-1);
else
    showlevels = round(sort(levels) * maxdens);
end

% Set up evenly spaced colors.
% if (exist('color'))
%     cmap = color;
% else
    cmap = colormap;
% end
cind = floor(linspace(1,size(cmap,1),length(showlevels)));

% Construct index vectors from the original ranges of the data.
x_inds = linspace(min1,max1,bins);
y_inds = linspace(min2,max2,bins);
z_inds = linspace(min3,max3,bins);

% Prepare figure and prettify.
%figreset;                    % too easy to get muddled if properties aren't defaults
set(gcf, 'Colormap', cmap);  % figreset changes this, so we reset as a courtesy.
axhand = gca;
grid on;                     % \ these two make it easier to follow when rotating the
set(axhand,'box','on');         % /         volume
daspect([1,1,1]);
view(240,20);                % no reason, just a generically good angle to start from

% Next we actually draw the contours.
for contour = 1:length(showlevels)
    fv = isosurface(x_inds,y_inds,z_inds,smoothdata,showlevels(contour));
    fv = reducepatch(fv, patch_reduction);
    p = patch(fv);
    
    % Color is spread over the colormap.
    set(p, 'FaceColor', cmap(cind(contour),:), 'EdgeColor', 'none');
    % The alpha is set to the reciprocal of the number of contours.
    set(p, 'AlphaDataMapping', 'none', 'FaceAlpha', 1/(length(showlevels)));

    % keep the axes looking the same in case 'progress' is turned on.
    axis tight;
    set(gca,'CameraViewAngleMode','manual', 'CameraViewAngle', 10);
    drawnow;
end

% Prepare the axes for rotating and add some lighting 
axis tight;
set(axhand,'CameraViewAngleMode','manual', 'CameraViewAngle', 10);
material dull;
lighting phong;
lightangle(90,60);
lightangle(270,10);

% Finally, color the tick marks.
set(axhand, 'Xcolor', [0.6 0 0]);
set(axhand, 'Ycolor', [0 0.6 0]);
set(axhand, 'Zcolor', [0 0 0.6]);
