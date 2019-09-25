% SBM Matlab utilities  
% Last updated: 09-15-2005
%
%
% DATATOOLS        [Matrix manipulation, signal processing and statistics]
%   cmtm              - Coherence function estimate via a multitaper method 
%   cmplx2rgb         - Convert complex data to RGB data via HSB space
%   convnsep          - Separable N-dimensional convolution
%   embed             - Lag embedding for scalar time series
%   gaussianvectors   - Matrix with multivariate Gaussian vectors as rows
%   gausskernel       - N-dimensional discretized Gaussian kernels 
%   hist3d            - 3-D density histogram
%   histnd            - N-D histogram
%   histxt            - Column-by-column histogram
%   histxy            - 2-D density histogram
%   knn               - K-Nearest Neighbors (naive algorithm)
%   leading_edges     - Finds 0->nonzero transitions
%   minmax            - Simultaneously find min/max elements of N-D arrays
%   mindist           - Indices of min distances in a pairwise dist mtx
%   padmatrix         - Add rows/columns to a matrix
%   pairdist          - Compute a Euclidean distance matrix
%   pcasvd            - Principal Components Analysis
%   pxcorr            - Cross-correlation estimate for point process data
%   rasterize         - Convert a point process to a binary time series
%   rescale           - Rescale the range of a data set
%   smooth3f          - Smooth 3D data (fast algorithm)
%   trig_events       - Extracts events from a time series
%
% DATATYPES        [Matlab data type manipulation]
%   isvector          - True for 1-D arrays
%   structcat         - Concatenates two structures field-by-field
%   structindex       - Indexes each field in a structure of arrays
%   structpack        - Copies workspace variables into structure field
%   structunpack      - Copy structure fields to workspace variables
%
% GRAPHICSTYLES    [Stylistic elements for visual consistency]
%   Cdgy              - Dark grey color vector
%   Clgy              - Light gray color vector
%   Cmar              - Maroorn color vector
%   Cmbl              - Medium blue color vector
%   Cmgr              - Medium green color vector
%   fancy             - Test OpenGL graphics
%   jetm              - Muted JET colormap
%   maplecolor        - Truecolor mimicking default Maple color scheme
%   nice              - Generic visual style settings
%   top               - Set camera to make a surface look like an image
%   xray              - Inverted grayscale colormap
%
% HGTOOLS          [Handle graphics UI menus and tools]
%   HGcolorshift      - Interactive Colorbar
%   HGlogmenu         - Toggles logarithmic scaling
%   linelabel         - Identify indices of lines plotted in a series
%
% MATLABTOOLS      [General utilities for Matlab programming]
%   busyfigure        - Disable/restore a UI figure
%   datetime2ind      - Convert date/time strings to serial time indices
%   printf            - Display formatted text
%   progressBar       - Improved progress indicator
%   ps2pdf            - Converts PostScript to Adobe PDF.
%   report_addpage    - Create multipage PostScript from Matlab figures
%   strmatchre        - Select strings matching regular expressions
%
% PLOTTYPES        [Graphical display of data]
%   anglehist         - Magnitude-weighted histogram of polar data
%   errorarea         - Line with confidence regions
%   linec             - Line with varying color
%   movie_player      - Displays 3-D matrix as a movie
%   mplot             - Efficient plotting of large numbers of lines
