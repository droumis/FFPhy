function [out, distances] = simulateRateModel(startdata, enddata, model)
%out = simulateRateDecline(startdata, model)
%
%MODEL - 'linear' or 'sigmoidal'

%values = values(find(~isnan(values)));
enddist = flipud(sortrows(enddata,1));
numvalues = size(enddist,1);
enddist(:,2) = 1;
enddist(:,2) = cumsum(enddist(:,2))/numvalues;
enddist = flipud(enddist);

binstep = .01;
Xi = [min(enddist(:,1)):binstep:max(enddist(:,1))]';
endcumsum = enddist(lookup(Xi,enddist(:,1)),2);

if isequal(model,'linear')
    steps = [0:.01:2];
    distances = [];
    for i = 1:length(steps)
        simdata = startdata*steps(i);
        simdist = flipud(sortrows(simdata,1));
        maxsimdist = max(simdist);
        numvalues = size(simdist,1);
        simdist(:,2) = 1;
        simdist(:,2) = cumsum(simdist(:,2))/numvalues;
        simdist = flipud(simdist);
        simcumsum = simdist(lookup(Xi,simdist(:,1)),2);
        simcumsum(find(Xi > maxsimdist)) = 0; 
        distances(i) = mean(abs(simcumsum-endcumsum));
    end
    [mindist, minind] = min(distances);
    
    out = [startdata*steps(lookup(.25,steps)) startdata*steps(lookup(.5,steps)) startdata*steps(lookup(1,steps)) startdata*steps(minind)];
    
elseif isequal(model,'sigmoid')
    slopesteps = [.01:.01:2];
    shift = [0:40];
    distances = [];
    for i = 1:length(slopesteps)
        for j = 1:length(shift)        
            scaling = 1./(1+exp(-slopesteps(i)*(Xi-shift(j))));
            simdata = startdata.* scaling(lookup(startdata,Xi));
            simdist = flipud(sortrows(simdata,1));
            maxsimdist = max(simdist);
            numvalues = size(simdist,1);
            simdist(:,2) = 1;
            simdist(:,2) = cumsum(simdist(:,2))/numvalues;
            simdist = flipud(simdist);
            simcumsum = simdist(lookup(Xi,simdist(:,1)),2);
            simcumsum(find(Xi > maxsimdist)) = 0;
            distances(i,j) = mean(abs(simcumsum-endcumsum));
        end
    end
    [mindist, minind] = min(distances(:));
    [indi,indj] = ind2sub([length(slopesteps) length(shift)],minind);
    optparam = [slopesteps(indi) shift(indj)];
    %scaling = 1./(1+exp(-slopesteps(indi)*(Xi-shift(indj))));
    scaling = 1./(1+exp(-.3*(Xi-15)));
    out = startdata.* scaling(lookup(startdata,Xi)); 
end
