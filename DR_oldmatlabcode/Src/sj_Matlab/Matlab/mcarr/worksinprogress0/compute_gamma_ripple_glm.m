function stats = compute_gamma_ripple_glm(directoryname,fileprefix,days,varargin)

%
%Goes through each day and computes the glm model for how predictive ca1
%and ca3 gamma power is to the presence of a ripple

%directoryname - example '/data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers 
%

% check to see if the directory has a trailing '/'
if (directoryname(end) ~= '/')
    warning('directoryname should end with a ''/'', appending one and continuing');
    directoryname(end+1) = '/';
end

x = -10:0.5:20;

stats.ca1_beta = [];
stats.ca1_pvalue = [];
stats.ca1_prediction = [];

stats.ca3_beta = [];
stats.ca3_pvalue = [];
stats.ca3_prediction = [];

%Normalize by the mean and std of the first session
mean_ca1 = []; mean_ca3 = []; std_ca1 = []; std_ca3 = [];
days = days(:)';
for day = days
    %load the ripples, tetinfo, and cellinfo files
    tetinfo = loaddatastruct(directoryname,fileprefix,'tetinfo');
    cellinfo = loaddatastruct(directoryname,fileprefix,'cellinfo');
    ripples = loaddatastruct(directoryname,fileprefix,'ripples');
    for epoch = 1:min(7,length(ripples{day}))
        if ~isempty(ripples{day}{epoch})

        % Determine which tetrodes are in CA1 and which are in CA3
        ca1tetrodes = evaluatefilter(tetinfo{day}{epoch},'isequal($area,''CA1'') & $numcells>1');
        ca3tetrodes = evaluatefilter(tetinfo{day}{epoch},'isequal($area,''CA3'') & $numcells>1');

        % Go through each tetrode and compute the z-score CA3 and CA1 lowgamma envelope
        ca1gam = []; ca3gam = [];
        for t = ca1tetrodes'
            loadfile = sprintf('%s/EEGnonreference/%slowgamma%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,t);
            load(loadfile)
            if isempty(ca1gam)
                ca1time = geteegtimes(lowgamma{day}{epoch}{t});
                ca1gam = double(lowgamma{day}{epoch}{t}.data(:,3));
                count = 1;
            else
                if size(lowgamma{day}{epoch}{t}.data(:,3)) == size(ca1gam)
                    ca1gam = ca1gam + double(lowgamma{day}{epoch}{t}.data(:,3));
                    count = count+1;
                elseif size(lowgamma{day}{epoch}{t}.data(:,3)) >= size(ca1gam)
                    tmptime = geteegtimes(lowgamma{day}{epoch}{t});
                    ind = lookup(ca1time,tmptime); clear tmptime
                    ca1gam = ca1gam + double(lowgamma{day}{epoch}{t}.data(ind,3));
                    count = count+1;
                else
                    tmptime = geteegtimes(lowgamma{day}{epoch}{t});
                    ind = lookup(tmptime,ca1time);
                    ca1time = tmptime; clear tmptime
                    ca1gam = ca1gam(ind) + double(lowgamma{day}{epoch}{t}.data(:,3));
                    count = count+1;
                end
            end
        end
        ca1gam = ca1gam./count;
        if isempty(mean_ca1)
            mean_ca1 = mean(ca1gam);
            std_ca1 = std(ca1gam);
        end
        ca1gam = (ca1gam - mean_ca1)./std_ca1;

        for t = ca3tetrodes'
            loadfile = sprintf('%s/EEGnonreference/%slowgamma%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,t);
            load(loadfile)
            if isempty(ca3gam)
                ca3time = geteegtimes(lowgamma{day}{epoch}{t});
                ca3gam = double(lowgamma{day}{epoch}{t}.data(:,3));
                count = 1;
            else
                if size(lowgamma{day}{epoch}{t}.data(:,3)) == size(ca3gam)
                    ca3gam = ca3gam + double(lowgamma{day}{epoch}{t}.data(:,3));
                    count = count+1;
                elseif size(lowgamma{day}{epoch}{t}.data(:,3)) >= size(ca3gam)
                    tmptime = geteegtimes(lowgamma{day}{epoch}{t});
                    ind = lookup(ca3time,tmptime); clear tmptime
                    ca3gam = ca3gam + double(lowgamma{day}{epoch}{t}.data(ind,3));
                    count = count+1;
                else
                    tmptime = geteegtimes(lowgamma{day}{epoch}{t});
                    ind = lookup(tmptime,ca3time);
                    ca3time = tmptime; clear tmptime
                    ca3gam = ca3gam(ind) + double(lowgamma{day}{epoch}{t}.data(:,3));
                    count = count+1;
                end
            end
        end
        ca3gam = ca3gam./count;
        if isempty(mean_ca3)
            mean_ca3 = mean(ca3gam);
            std_ca3 = std(ca3gam);
        end
        ca3gam = (ca3gam - mean_ca3)./std_ca3;
        clear lowgamma

        riptimes = getripples([day epoch], ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1''))','minstd',3);
        if isequal(ca1time,ca3time)
            time = ca1time;
            clear ca1time ca3time
        elseif length(ca1time)>length(ca3time)
            ind = lookup(ca3time,ca1time);
            time = ca3time;
            ca1gam = ca1gam(ind);
            clear ca1time ca3time ind
        else
           ind = lookup(ca3time,ca1time);
            time = ca1time;
            ca3gam = ca3gam(ind);
            clear ca1time ca3time ind
        end

        times = zeros(size(time));
        for i = 1:size(riptimes,1)
            ind = lookup(riptimes(i,[1 2]),time);
            times(ind(1):ind(2)) = 1;
        end
        
        %Subdivide the data into 200 ms windows
        bin = round(0.2/median(diff(time)));
        subs = lookup(time,time(1:bin:end));
        ca1gam = accumarray(subs,ca1gam,[max(subs) 1],@(x) nanmax(x));
        ca3gam = accumarray(subs,ca3gam,[max(subs) 1],@(x) nanmax(x));
        times = accumarray(subs,times,[max(subs) 1],@(x) nanmax(x));
       
        %Compute GLM model to determine significance
        [beta,dev,tmpstats] = glmfit(ca1gam,times,'binomial');

        %Make histogram of data to show where gamma power predicts ripples
        subs = lookup(ca1gam,x);
        pred1 = accumarray(subs,times,[length(x) 1],@(x) mean(x), NaN);
        
        stats.ca1_beta = [stats.ca1_beta beta];
        stats.ca1_pvalue = [stats.ca1_pvalue tmpstats.p];
        stats.ca1_prediction = [stats.ca1_prediction pred1];
        
        %Compute GLM model to determine significance
        [beta,dev,tmpstats] = glmfit(ca3gam,times,'binomial');

        %Make histogram of data to show where gamma power predicts ripples
        subs = lookup(ca3gam,x);
        pred3 = accumarray(subs,times,[length(x) 1],@(x) mean(x), NaN);
        
        stats.ca3_beta = [stats.ca3_beta beta];
        stats.ca3_pvalue = [stats.ca3_pvalue tmpstats.p];
        stats.ca3_prediction = [stats.ca3_prediction pred3];
        end
    end
end

            