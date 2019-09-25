function out = calcprobdist(values,xvalues)

%plot probability distribution plots
%-------------------------------------------

%binstep = .01;
smoothwindow = round(length(xvalues)/30);
values = values(find(~isnan(values)));
values = flipud(sortrows(values,1));
numvalues = size(values,1);
values(:,2) = 1;
values(:,2) = cumsum(values(:,2))/numvalues;
[Xval, Xindex] = unique(values(:,1));
Yval = values(Xindex,2);

%Xi = [min(Xval):binstep:max(Xval)]';
Xi = xvalues;
Yi = Yval(lookup(Xi,Xval));

%Yi = interp1(Xval,Yval,Xi);

Yslope = diff(Yi);
Yslope = [Yslope(1);Yslope(:)];

Xout = Xi(2:end);
Yout = -smoothvect(Yslope,gaussian(smoothwindow,smoothwindow*10));

out = [Yout(:)/(sum(Yout))]';

    
    




