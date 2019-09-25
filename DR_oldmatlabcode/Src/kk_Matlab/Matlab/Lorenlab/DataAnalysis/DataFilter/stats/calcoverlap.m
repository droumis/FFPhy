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
    if length(trajdata1)>=traj & length(trajdata2)>=traj
        if ~isempty(trajdata1{traj}) & ~isempty(trajdata2{traj})
            peaks = [peaks max([max(trajdata1{traj}(:,5)) max(trajdata2{traj}(:,5))])];
        end
    end
end
peak = max(peaks);
overlap = [];
if (peak >= thresh)
    for traj = 1:length(trajdata1)
        if length(trajdata1)>=traj & length(trajdata2)>=traj
            if ~isempty(trajdata1{traj}) & ~isempty(trajdata2{traj})
                trajlengths = [length(trajdata1{traj}(:,5)) length(trajdata2{traj}(:,5))];

                trajlength = min(trajlengths);
                if normalize==1
                    ratediff = abs((trajdata1{traj}(1:trajlength,5)/max(trajdata1{traj}(1:trajlength,5))) - (trajdata2{traj}(1:trajlength,5)/max(trajdata2{traj}(1:trajlength,5))));
                    totalrates = ((trajdata1{traj}(1:trajlength,5)/max(trajdata1{traj}(1:trajlength,5))) + (trajdata2{traj}(1:trajlength,5)/max(trajdata2{traj}(1:trajlength,5))));
                elseif normalize ==0
                    ratediff = abs((trajdata1{traj}(1:trajlength,5)) - (trajdata2{traj}(1:trajlength,5)));
                    totalrates = ((trajdata1{traj}(1:trajlength,5)) + (trajdata2{traj}(1:trajlength,5)));
                end
                ratediff(find(isnan(ratediff))) = 0;
                totalrates(find(isnan(totalrates))) = 0;
                if (sum(totalrates) > 0)
                    overlap = [overlap (sum(totalrates)-sum(ratediff))/(sum(totalrates))]; %overlap;
                end
            end
        end

    end
end
if ~isempty(overlap)
    overlap = mean(overlap);
else %no firing in either traj or one traj isempty
    overlap = NaN;
end