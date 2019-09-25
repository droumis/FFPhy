function H= behaventropy(x, y, xedges, yedges, times, excludeperiods)
% takes xy position and bin parameters and determines behavioral entropy:
% the entropy of the transitions between spaitial positions.  
% H= sum over i & j( -prob of transition from j to i * log2(prob of
% transition from j to i) )
% Higher entropy = more variable behavior.
%
% see Jackson, Johnson, & Redish J Neuroscience 2006
%
% INPUTS: 
%   x: vector of x positions, ie 2nd column of pos structure
%   y: vector of y positions, ie 3nd column of pos
%   xedges: bin edges along x axis
%   yedges: bin edges along y axis
%   times: time of each position, row corresponds to row of x & y, ie 1st
%   column of pos
%   excludeperiods: if empty no times will be excluded
%
% OUTPUT: behavioral entropuy
%
% AS 4/15/09


%turn into periods-type matrix with start col and end col, have to add to
%start and end col so that all values that are eual to a bin edge fall into
%a bin
xedges = [xedges(1:end-1)+.0001; xedges(2:end)+.00001]';
yedges = [yedges(1:end-1)+.0001; yedges(2:end)+.00001]';
xbinind = getbinindex(x,xedges);
ybinind = getbinindex(y,yedges);
if any(xbinind ==0) || any(ybinind ==0)
    Error('Position falls outside of bins')
end

%turn xy index into other index
xy = reshape([1: size(yedges,1) * size(xedges, 1)], size(yedges,1), size(xedges, 1)); %not used but useful to visualize bins
xybin = sub2ind([size(yedges,1) size(xedges,1)], ybinind(ybinind~=0 & xbinind~=0), xbinind(xbinind~=0 & ybinind~=0) );
% one epoch starts with 2 points with negative x --> above line excludes those
% points

%find where bin change
binchind = (find(diff(xybin)));
binch = [xybin(binchind)  xybin(binchind+1) times(binchind) ]; %transitions: current bin, bin goes to

if ~isempty(excludeperiods)
goodbinch = binch(~isExcluded( binch(:,3), excludeperiods),:); %what rows should be included of binch, ie fall into includetimes
else
    goodbinch = binch;
end
%compuyte transition probabilities for each transition
st = unique(goodbinch(:,1)); %current bins/start bins for each bin change
allp=[];
for a= 1:length(st)
    trans = goodbinch(goodbinch(:,1)==st(a),:);
    table = tabulate(trans(:,2));
    probs = table(:,3)*.01;
    probs = probs(probs>0); %if prob = 0, exclude from anal
    allp = [allp;probs];
end
H = sum(-allp .* log2(allp)); %H = behaventropy(x, y, xedges, yedges);


