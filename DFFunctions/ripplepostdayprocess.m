function outliers = ripplepostdayprocess( animdirectory,fileprefix, days, peakthresh, varargin)
% RIPPLEPOSTDAYPROCESS(directoryname,fileprefix,days, options)
% options -
%   'exclude' , 1 or 0, default 0
%           1: exclude outlier ripples
% PEAKTHRESH peaks over peakthresh will be considered outlier ripples
% OUTLIERS list of outlier ripples per epoch: [day epoch tetrode #outliers
% ripples.thresh]
%
% run after running extractripples to identify and exclude outlier ripples with peak >
% peakthresh

% assign the options
exclude = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'exclude'
            exclude = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end


outliers = [];
for i = 1:length(days)
    d = days(i);
    dsz = '';
    if d <10
        dsz = '0';
    end
    eval(['load ', animdirectory, fileprefix, 'ripples', dsz, num2str(d)]);
    for e = 1:length(ripples{d}) %for all epochs in ripples{d}
        for t = 1:length(ripples{d}{e}) %for all tetrodes
            if ~isempty(ripples{d}{e}{t})
                out = findrippleoutliers([d e t], [], ripples, peakthresh, 'appendindex', 1);
                %find outlier ripples and their index
                if out(:,4) > 0  %if there were some outlier ripples add them to the output
                    outliers = [outliers; out];

                    if (exclude) %exclude outlier ripples
                        rip = ripples{d}{e}{t};
                        include = find(rip.peak <= peakthresh);
                        n = fieldnames(rip);
                        for i =  1: 10 %for each field, select only included ripples
                            eval(['ripples{d}{e}{t}.' n{i} ' = rip.' n{i} '(include)'])
                        end
                    end

                end
            end
        end
    end
    if (exclude)
    eval(['save ',animdirectory,fileprefix,'ripples',dsz,num2str(d),' ripples']);
    %save new ripple structure
    end
end
end