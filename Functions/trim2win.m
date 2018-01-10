

function out = trim2win(inData, srate, plotwin)

midpoint = (size(inData,2)-1)/2; %get middle index of window
out = inData(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot wind
end