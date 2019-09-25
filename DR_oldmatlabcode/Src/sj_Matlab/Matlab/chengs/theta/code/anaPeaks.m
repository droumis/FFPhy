function [ix, lx, hy, wx]= anaPeaks(y,x,widthRatio, sideId)
%function [hy, lx, wx]= anaPeaks(y,x)
%  Find maxima in y and return 
%  ix: indices of peaks
%  lx: location of peaks
%  hy: height of peaks
%  wx: width at widthRatio of max
%
%  Really should do local quadratic fit and determine x of max. of fitted function.

if nargin<3; widthRatio= 0.5; end
if nargin<4 | strcmp(sideId, 'both')
    side= 2;  % both sides
else
    if strcmp(sideId, 'left')
        side= 0;
    elseif strcmp(sideId, 'right')
        side= 1;
    else 
        error('unknown sideId');
    end
end


n= length(y);
if size(y,1)==1
    yl= [inf, y(1:n-1)];
    yr= [y(2:n), inf];
else
    yl= [inf; y(1:n-1)];
    yr= [y(2:n); inf];
end
ix= find(yl < y & yr < y);

lx= x(ix);
hy= y(ix);
wx= zeros(size(ix));

for iix=1:length(ix)
    i= ix(iix);
    cut= widthRatio*y(i);
    if side==0 | side==2
        j= i;
        while j>=1 & y(j)>cut; 
%            if(y(j)>y(i)); j=1; break; end
            j= j-1; 
        end
        if j>0; 
            xlo= interp1(y(j:j+1),x(j:j+1),cut);
        else 
            xlo= nan; 
        end
    end

    if side==1 | side==2
        j= i;
        while j<= n & y(j)>cut; 
%            if(y(j)>y(i)); j=n; break; end
            j= j+1; 
        end
        if j<=n; 
            xhi= interp1(y(j-1:j),x(j-1:j),cut);
        else 
            xhi= nan; 
        end
    end
%    keyboard

    switch side
    case 0
        wx(iix)= x(i)-xlo;
    case 1
        wx(iix)= xhi-x(i);
    case 2
        wx(iix)= xhi-xlo;
    end
end
