function out = tally(x);
%TALLY	TALLY(X) calculates the frequencies of the distinct levels of X.
%	Output is a matrix with two rows, the first containing the
%	levels of X and the second the corresponding counts.

%	GKS  27 Dec 1992

s=sort(x(:))';
[m n]=size(s);
i=find( [ 1 diff(s) ] > 0 );
out=[ s(i); diff([i n+1]) ];
