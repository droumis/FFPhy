function smoothhist2D(X,lambda,nbins,outliercutoff,plottype)
% SMOOTHHIST2D Plot a smoothed histogram of bivariate data.
%   SMOOTHHIST2D(X,LAMBDA,NBINS) plots a smoothed histogram of the bivariate
%   data in the N-by-2 matrix X.  Rows of X correspond to observations.  The
%   first column of X corresponds to the horizontal axis of the figure, the
%   second to the vertical. LAMBDA is a positive scalar smoothing parameter;
%   higher values lead to more smoothing, values close to zero lead to a plot
%   that is essentially just the raw data.  NBINS is a two-element vector
%   that determines the number of histogram bins in the horizontal and
%   vertical directions.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF) plots outliers in the data as points
%   overlaid on the smoothed histogram.  Outliers are defined as points in
%   regions where the smoothed density is less than (100*CUTOFF)% of the
%   maximum density.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,[],'surf') plots a smoothed histogram as a
%   surface plot.  SMOOTHHIST2D ignores the CUTOFF input in this case, and
%   the surface plot does not include outliers.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF,'image') plots the histogram as an
%   image plot, the default.
%
%   Example:
%       X = [mvnrnd([0 5], [3 0; 0 3], 2000);
%            mvnrnd([0 8], [1 0; 0 5], 2000);
%            mvnrnd([3 5], [5 0; 0 1], 2000)];
%       smoothhist2D(X,5,[100, 100],.05);
%       smoothhist2D(X,5,[100, 100],[],'surf');
%
%   Reference:
%      Eilers, P.H.C. and Goeman, J.J (2004) "Enhancing scaterplots with
%      smoothed densities", Bioinformatics 20(5):623-628.

%   Copyright 2009 The MathWorks, Inc.
%   Revision: 1.0  Date: 2006/12/12
%
%   Requires MATLAB� R14.


% If either the number or arguments is less than 4 or the outliers field 
% is empty, the set outlier value to 0.5 
if nargin < 4 || isempty(outliercutoff), outliercutoff = .05; end
% If the number of arguments is less than 5 plot image
if nargin < 5, plottype = 'image'; end

% --------------------------------------
% First part transforms data into a reference frame with a specified number
% of frames (nbins)

% Convert data to log

% Generates boundries of the frames
% Finds minimum number in the elements of the input matrix called X  
minx = min(X,[],1);
% Finds maximum number in the elements of the input matrix called X 
maxx = max(X,[],1);
rminx = ceil(minx);
rmaxx = floor(maxx);

% X axis
% Divides up the the range of input data into bins specified by nbins
edges1 = linspace(rminx(1), rmaxx(1), (rmaxx(1)-rminx(1))/nbins(1)+1);
% Creates the midway point between values specified in edges1
ctrs1 = edges1(1:end-1) + .5*diff(edges1);
% Creates vector from data with divisions of number of bins and infinity at
% the end
edges1 = [-Inf edges1(2:end-1) Inf];

% Y axis
% Does the same for the other dimension of input data
edges2 = linspace(rminx(2), rmaxx(2), (rmaxx(2)-rminx(2))/nbins(2)+1);
ctrs2 = edges2(1:end-1) + .5*diff(edges2);
edges2 = [-Inf edges2(2:end-1) Inf];

% Transforms data into reference grid dimensions
% Makes a matrix with the same size as input data
[n,p] = size(X);

% Makes a matrix with n rows with 2 columns
bin = zeros(n,2);
% Reverse the columns of H to put the first column of X along the
% horizontal axis, the second along the vertical.
% Generates a matrix that has a cumulative count of how many values there
% are inbetween the intervals defined in edges1
[dum,bin(:,2)] = histc(X(:,1),edges1);
[dum,bin(:,1)] = histc(X(:,2),edges2);

% ----------------------------------------

% Sums up values belong to each pixel in the new reference frame
% Perform classification on matrix bin to fit into a array the size
% specified by nbins
xdim = (rmaxx(2)-rminx(2))/nbins(2)+1;
ydim = (rmaxx(1)-rminx(1))/nbins(1)+1;
E = accumarray(bin,1,[xdim ydim]);
H = E ./ n;

% Eiler's 1D smooth, twice
G = smooth1D(H,lambda);
F = smooth1D(G',lambda)';
% % An alternative, using filter2.  However, lambda means totally different
% % things in this case: for smooth1D, it is a smoothness penalty parameter,
% % while for filter2D, it is a window halfwidth
% F = filter2D(H,lambda);

%expresses results of F as a fraction of its maximum value
relF = F./max(F(:));
if outliercutoff > 0
    outliers = (relF(nbins(2)*(bin(:,2)-1)+bin(:,1)) < outliercutoff);
end

% Sets color scale with number of colors
nc = 256;
% Chooses color scheme see available matlab schemes
% colormap(hot(nc));
colormap(jet(nc));




% Chooses the type of plot made

switch plottype
case 'surf'
    surf(ctrs1,ctrs2,F,'edgealpha',0);
case 'image'
    % Makes image with dimensions ctrs1 and ctrs2 and each pixel has value
    % of relF spanning the scale defined by nc
    imagetitl='Run 1 positions'
    im=image(ctrs1,ctrs2,nc.*relF); title(imagetitl);
    % sets the image to lock original proportions
    axis image;
    
    grid on;
    
    set(gca,'YTick',[min(ctrs2):10:max(ctrs2)]);
    
    xlabel('cm');
    ylabel('cm');
    
   
    % Following part does not work yet
    % Selects lowest cutoff
    %lowF = 0.1*max(F(:));
    % Selects high cutoff
    %highF = 0.5*max(F(:));
    % Creates colormap between values
    %myscale=[lowF highF];
    %caxis manual;
    %caxis(myscale);

    
    % Adds color legend, tick names
    % Number of ticks between min and max values
    scaleint=[5];
    % A list of numbers corresponding to the intervals of the colorbar
    datascale=[0:nc/scaleint:nc];
    % A list of numbers corresponding to the names of the intervals
    tickname = round([min(E(:)):max(E(:))/scaleint:max(E(:))]);
    % Make the colorbar, specifying the location and the tick names
    c=colorbar('location','EastOutside','YTickLabel', tickname);
    % Adds text 'Legend' to color legend
    ylabel(c,'Events');
    % Changes the tick intervals of color legend at intervals defined by scalteint
    set(c, 'YTick',datascale);
    
    
    
    % image(ctrs1,ctrs2);
    hold on
    % plot the outliers in grey
    if outliercutoff > 0
        plot(X(outliers,1),X(outliers,2),'.','MarkerEdgeColor',[.8 .8 .8]);
        
     
    end
%     % plot a subsample of the data
%     Xsample = X(randsample(n,n/10),:);
%     plot(Xsample(:,1),Xsample(:,2),'bo');
    hold off
end

%-----------------------------------------------------------------------------
function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
% This is a better solution, but takes a bit longer for n and m large
% opts.RECT = true;
% D1 = [diff(E,1); zeros(1,n)];
% D2 = [diff(D1,1); zeros(1,n)];
% Z = linsolve([E; 2.*sqrt(lambda).*D1; lambda.*D2],[Y; zeros(2*m,n)],opts);


%-----------------------------------------------------------------------------
function Z = filter2D(Y,bw)
z = -1:(1/bw):1;
k = .75 * (1 - z.^2); % epanechnikov-like weights
k = k ./ sum(k);
Z = filter2(k'*k,Y);
