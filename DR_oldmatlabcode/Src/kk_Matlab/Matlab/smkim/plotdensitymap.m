function h = plotdensitymap(d,occ,dims,occthresh)
%PLOTDENSITYMAP Plots a density estimate in 1 or 2 dimensions.
%   PLOTDENSITYMAP(D,OCC,DIMS,OCCTHRESH) plots a 1- or 2- 
%   dimensional density estimate D as a colormap. Alpha  
%   opacity is set by the occupancy OCC and OCCTHRESH. 
%
%   DIMS is a struct array of length N<=2, with the 
%   following fields:
%       NAME is a human-readable descriptive string 
%       (e.g.'distance (cm)').
%       TYPE is 'linear' or 'circular'
%       GRID is a vector of evenly-spaced grid points at 
%       which the kernel density is to be estimated; if 
%       TYPE is 'circular', then GRID(1) coincides with 
%       GRID(end)
%   
%   OCC and D must be matching vectors or matrices of size 
%   [CELLFUN(@LENGTH,{DIMS(:).GRID})-1, filled with 
%   nonnegative real values.
%
%   See also DENSITYMAP.
%
% Written by smk, 7 November 2008
%

n = ndims(d);
if n > 2
    error('D and OCC may must be 1- or 2-dimensional');
end
if any(size(d) ~= size(occ))
    error('D and OCC must be the same size');
end
if sireal(d) || ~isreal(occ) || any(d(:) < 0) || any(occ(:) < 0)
    error('D and OCC must contain only non-negative real values');
end
if ~isstruct(dims) || ((length(dims) ~= 1) && (length(dims) ~= 2))
    error('DIMS must be a struct array of length no greater than 2');
end
if length(dims) ~= n
    error('dimensionality of D and OCC must equal length of DIMS');
end
for i = 1:n
    if any(size(d) ~= cellfun(@length,{dims(:).grid})-1)
        error('size of D does not match the grid defined in DIMS');
    end
    if all(diff(diff(dims(i).grid)) ~= 0)
        error('grid spacings in DIMS must be uniform monotonic');
    end
end

% set alpha values to censor the density plot according to
% occupancy density
a = occ;
a(occ < 0) = 0;
a(occ > occthresh) = 1;

if n > 2
    error('can not plot a density map in higher than 2 dimensions!')
elseif n==2
    % use IMAGE to plot pixels
    xgrid=0.5*(dims(1).grid(1:end-1)+dims(1).grid(2:end));
    ygrid=0.5*(dims(2).grid(1:end-1)+dims(2).grid(2:end));
    h = image(xgrid,ygrid,d','AlphaData',a', ...
      'CDataMapping','scaled');
    view([0 90]);
    set(gca,'Box','on','DataAspectRatio',[1 1 1], ...
      'CLim',[min(d(a>=1)) max(d(a>=1))],'YDir','normal');
elseif n==1
    % plot as a LINE object
    xgrid=0.5*(dims(1).grid(1:end-1)+dims(1).grid(2:end));
    ygrid = d;
    ygrid(a < 1) = NaN;
    disp(nnz(isnan(ygrid)));
    h = line('XData',xgrid,'YData',ygrid,'LineStyle','-');
    set(gca,'XLim',[dims(1).grid(1) dims(1).grid(end)]);
else
    error('DIMS is not valid');
end
