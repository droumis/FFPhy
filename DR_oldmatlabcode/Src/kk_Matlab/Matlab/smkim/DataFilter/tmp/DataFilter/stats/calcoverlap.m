function overlap = calcoverlap(trajdata1,trajdata2, varargin)


normalize = 0;
thresh = 0;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'normalize'
                normalize = varargin{option+1};
            case 'thresh'
                thresh = varargin{option+1};    
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

peaks = [];
for traj = 1:length(trajdata1)
    peaks = [peaks max([max(trajdata1{traj}(:,5)) max(trajdata2{traj}(:,5))])];
end
peak = max(peaks);
overlap = [];
ratediff = [];
totalrates = [];

for traj = 1:length(trajdata1)
    trajlengths = [length(trajdata1{traj}(:,5)) length(trajdata2{traj}(:,5))];

    trajlength = min(trajlengths);
    %ratediff = abs((trajdata1{traj}(1:trajlength,5)/max(trajdata1{traj}(1:trajlength,5))) - (trajdata2{traj}(1:trajlength,5)/max(trajdata2{traj}(1:trajlength,5))));
    %totalrates = ((trajdata1{traj}(1:trajlength,5)/max(trajdata1{traj}(1:trajlength,5))) + (trajdata2{traj}(1:trajlength,5)/max(trajdata2{traj}(1:trajlength,5))));

    ratediff = [ratediff; abs((trajdata1{traj}(1:trajlength,5)) - (trajdata2{traj}(1:trajlength,5)))];
    totalrates = [totalrates; ((trajdata1{traj}(1:trajlength,5)) + (trajdata2{traj}(1:trajlength,5)))];
end

ratediff(find(isnan(ratediff))) = 0;
totalrates(find(isnan(totalrates))) = 0;
if (max(totalrates) > thresh)
    overlap = (sum(totalrates)-sum(ratediff))/(sum(totalrates));
    %overlap = [overlap (sum(totalrates)-sum(ratediff))/(sum(totalrates))]; %overlap;
else
    overlap = []; %overlap;
end

