function stats = compute_synchrony_ripple_glm(directoryname,fileprefix,days,varargin)

%Goes through each day and computes the glm model for how predictive ca1
%and ca3 gamma synchrony is to the presence of a ripple

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

x = 0:0.1:1;

stats.beta = [];
stats.pvalue = [];
stats.prediction = [];

params = [];
params.Fs = 1500;
params.fpass = [1 350];
window =[0.2 0.2];


days = days(:)';
for day = days
    %load the ripples, tetinfo, and cellinfo files
    tetinfo = loaddatastruct(directoryname,fileprefix,'tetinfo');
    cellinfo = loaddatastruct(directoryname,fileprefix,'cellinfo');
    ripples = loaddatastruct(directoryname,fileprefix,'ripples');
    for epoch = 1:min(7,length(ripples{day}))
        if ~isempty(ripples{day}{epoch})

        % Determine which tetrodes are in CA1 and which are in CA3
        ca1tetrode = evaluatefilter(tetinfo{day}{epoch},'isequal($area,''CA1'') & $maxcell==1');
        ca3tetrode = evaluatefilter(tetinfo{day}{epoch},'isequal($area,''CA3'') & $maxcell==1');

            if ~isempty(ca1tetrode) && ~isempty(ca3tetrode)
                % Go through each tetrode pair and compute the coherence in 200 ms
                % windows
                loadfile = sprintf('%s/EEGnonreference/%seeg%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,ca1tetrode);
                load(loadfile)
                e1 = eeg{day}{epoch}{ca1tetrode};

                loadfile = sprintf('%s/EEGnonreference/%seeg%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,ca3tetrode);
                load(loadfile)

                e2 = eeg{day}{epoch}{ca3tetrode};
                e1start = e1.starttime;

                e1times = geteegtimes(e1);
                e2times = geteegtimes(e2);

                if length(e1times)>length(e2times)
                    temp = lookup(e2times,e1times);
                    e1 = e1.data(temp);
                    e2 = e2.data;
                    time = e2times;
                elseif length(e2times)>length(e1times)
                    temp = lookup(e1times,e2times);
                    e1 = e1.data;
                    e2 = e2.data(temp);
                    time = e1times;
                elseif length(e1times)==length(e2times)
                    e1 = e1.data;
                    e2 = e2.data;
                end

                [C,phi,S12,S1,S2,t,f]=cohgramc(e1,e2,window, params);
                %Figure out average gamma coherence
                gamind = lookup([20 50],f);
                gam = mean(C(:,gamind(1):gamind(2)),2);

                riptimes = getripples([day epoch], ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1''))','minstd',5);

                time = zeros(size(t));
                for i = 1:size(riptimes,1)
                    ind = lookup(riptimes(i,[1 2]),t+e1start);
                    time(ind(1):ind(2)) = 1;
                end

                %Compute GLM model to determine significance
                [beta,dev,tmpstats] = glmfit(gam,time','binomial');

                %Make histogram of data to show where gamma coherence predicts ripples
                subs = lookup(gam,x);
                pred1 = accumarray(subs,time,[length(x) 1],@(x) mean(x), NaN);

                stats.beta = [stats.beta beta];
                stats.pvalue = [stats.pvalue tmpstats.p];
                stats.prediction = [stats.prediction pred1];
            elseif isempty(stats.beta)
                stats.beta = [NaN; NaN];
                stats.pvalue = [NaN; NaN];
                stats.prediction = nan(size(x))';
            else
                stats.beta(:,end+1) = NaN;
                stats.pvalue(:,end+1) = NaN;
                stats.prediction(:,end+1) = NaN;
            end
        end
    end
end

            