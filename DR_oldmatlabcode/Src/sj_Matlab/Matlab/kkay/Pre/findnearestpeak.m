function index = findnearestpeak(data,startindex,searchwin)

%% returns index of data peak closest to startindex
%% uses first and second diff to find nearest local peak from startindex
%% searchwin is in samples, looks in window from -searchwin to +searchwin

datawindow=data((startindex-searchwin):(startindex+searchwin));    %% copies the small datawindow

rightdiff=((datawindow(1:(end-1))-datawindow(2:end))>0);
leftdiff=((datawindow-[0 datawindow(1:(end-1))])>0);
leftdiff=leftdiff(1:(end-1));

peak_indices=find(rightdiff & leftdiff);

dummy=0;
for i=peak_indices       %% iterate through all peaks to find the closest to startindex
    if abs(i-(searchwin+1)) < abs(dummy-(searchwin+1))
        dummy=i;
    end
end

index=startindex+dummy-(searchwin+1);

end