function [z,X]=MultNormal(mu,Sigma,borders)

% return multivariate normal PDF with mean mu and covariance Sigma between the 
% values specified by borders
%
% mu        row vector (with d entries)
% Sigma     d x d covariance matrix
% borders   [optional] (N+1) x d matrix, where each column gives the borders 
%           of the bin
%
% z         N^d array with probability densities

d=length(mu);
[d1 d2]=size(Sigma);
if d1 ~= d2
    error('Sigma must be square matrix');
end
if d1 ~= d | d2 ~= d
    error('Sigma''s dim must match mu''s');
end

% default borders
nBins=2;
nDev= 3;
if nargin <= 2
    dev=sqrt(diag(Sigma));
    for j=1:d
	N(j)=2*nBins*nDev+2;
	borders{j}= [-nDev-1/(2*nBins):1/nBins:nDev+1/(2*nBins)]'*dev(j) + ...
	    ones(N(j),1)* mu(j); 
    end
elseif length(borders)~=d
    errors('borders must have d elements');
end

for j=1:d
    X{j}= (borders{j}(1:end-1)+borders{j}(2:end))/2;
    N(j)= length(X{j});
end

Xpre=[];
Nd=1;
for j=1:d
    Xd=[];
    for k=1:N(j)
	Xd=[Xd; [Xpre, X{j}(k)*ones(Nd,1)]];
    end
    Nd=Nd*N(j);
    Xpre=Xd;
end
		     
Y=mvnpdf(Xd,mu,Sigma);
z=reshape(Y,N);
