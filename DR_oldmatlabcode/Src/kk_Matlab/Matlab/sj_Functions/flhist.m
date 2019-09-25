% [counts bincenters] = flhist(data, binedges)
% Calculate histogram.  Counts is the counts in each bin.  
% bincenter is the centers of each bin.
%
%  DATA is treated as a single array
%
%  binedges is a vector of bin boundaries (NOTE: this is different to the
%  mathworks HIST which takes bin centers)  If BINS is a scalar it is 
%  taken to be a bin count.
%
%  If BINS is specified as a count, or as a regularly spaced vector, a
%  rapid binning algorithm is used that depends on this and is linear
%  in ndata.  If BINS is irregularly spaced the algorithm used takes
%  time proportional to ndata*nbins.

