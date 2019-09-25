function xc= spikeCCG(t, M, coef, valid)
%function xc= spikeCCG(t, M, coef, valid)
%
% Use quick method to compute cross-correlations or autocorrelations

if nargin<3; coef= 1; end

% build aux variables
for i=1:length(t)
    if isempty(t{i}) | t{i}(end)-t{i}(1)<M; 
        warning('requested more lags than available data'); 
    end
    [su{i}, mi, mj]= unique(t{i});
    nsp= length(su{i});
    nu{i}= ones(1,nsp);
    for j=1:nsp; nu{i}(j)= sum(mj==j); end
end

if length(t)==1; 
    su{2}= su{1};
    nu{2}= nu{1};
end

% compute correlations
for m=-M:M
    [cv, cj1, cj2]= intersect(su{1}, su{2}+m);
    if nargin==4;
        ind= find(valid(cv) | valid(cv-m));
        cj1= cj1(ind);
        cj2= cj2(ind);
    end
    xc(m+M+1)= sum(nu{1}(cj1).*nu{2}(cj2));
end

% normalize
if coef
    for i=1:2
        if nargin==4;
            ind= find(valid(su{i}));
            n(i)= sum(nu{i}(ind).^2);
        else
            n(i)= sum(nu{i}.^2);
        end
    end
    if (n(1)*n(2))
        xc= xc/sqrt(n(1)*n(2));
    end
end
