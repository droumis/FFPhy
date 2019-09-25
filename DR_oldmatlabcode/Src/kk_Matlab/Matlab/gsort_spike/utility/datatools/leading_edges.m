function markers = leading_edges(data)
%LEADING_EDGES     Marks 0 -> NONZERO transitions in binary data.
%   EDGEMARKERS = LEADING_EDGES(BINARYDATA) takes an M x N matrix
%   BINARYDATA and returns a M x N matrix EDGEMARKERS, containing a '1' in
%   each location that corresponds to the first nonzero value in a
%   BINARYDATA column following one or more zeros above it in the same
%   column.  EDGEMARKERS is of type logical.
%
%   NOTE: If the first row of any BINARYDATA column contains a nonzero
%   value, it will be marked in EDGEMARKERS.

markers = [diff([repmat(logical(0), 1, size(data,2)); (data>0)], 1, 1) > 0];