function out = control_gammapower_spikemodulation(directoryname,fileprefix,days,varargin)

%
%Goes through each day and computes the gamma modulation of CA3 and CA1
%spikes for gamma power matched times with and without an SWR
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

out = cell(max(days),1);

%Normalize by the mean and std of the first session
mean_ca1 = []; std_ca1 = [];
days = days(:)';
for day = days
    %load the ripples, tetinfo, and cellinfo files
    tetinfo = loaddatastruct(directoryname,fileprefix,'tetinfo');
    cellinfo = loaddatastruct(directoryname,fileprefix,'cellinfo');
    ripples = loaddatastruct(directoryname,fileprefix,'ripples');
    spikes = loaddatastruct(directoryname,fileprefix,'spikes');
    
    out{day} = cell(min(7,length(ripples{day})),1);
    for epoch = 1:min(7,length(ripples{day}))
        
        %initialize output
        out{day}{epoch}.ripple_ca1_mod = [];
        out{day}{epoch}.ripple_ca3_mod = [];
        out{day}{epoch}.no_ripple_ca1_mod = [];
        out{day}{epoch}.no_ripple_ca3_mod = [];
        
        if ~isempty(ripples{day}{epoch})

            % Determine which tetrodes are in CA1 and which are in CA3
            ca1tetrodes = evaluatefilter(tetinfo{day}{epoch},'isequal($area,''CA1'')');
            ca3tetrodes = evaluatefilter(tetinfo{day}{epoch},'isequal($area,''CA3'') & $numcells>1');

            % Go through each tetrode and compute the z-score CA1 lowgamma envelope
            ca1gam = [];
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

            riptimes = getripples([day epoch], ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1''))','minstd',3);
            time = ca1time; clear ca1time

            times = zeros(size(time));
            for i = 1:size(riptimes,1)
                ind = lookup(riptimes(i,[1 2]),time);
                times(ind(1):ind(2)) = 1;
            end

            %Subdivide the data into 200 ms windows
            bin = round(0.2/median(diff(time)));
            subs = lookup(time,time(1:bin:end));
            tmpgam = zeros(max(subs),1); tmpind = zeros(max(subs),1);
            for i = 1:max(subs)
                [tmpgam(i) tmpind(i)] = max(ca1gam(subs==i));
                tmpind(i) = time(find(subs==i,1)-1+tmpind(i));
            end
            ca1gam = [tmpgam tmpind];
            times = accumarray(subs,times,[max(subs) 1],@(x) nanmax(x));

            %Make histogram of data to show where gamma power predicts ripples
            subs = lookup(ca1gam(:,1),x);
            pred1 = accumarray(subs,times,[length(x) 1],@(x) mean(x), NaN);

            %Go through and pick out 200 ms bins that correspond to 40-60%
            %predictability. Divide into rip+ and rip- times.
            if max(pred1 > 0.6)
                x_interp = -10:0.1:20;
                ca1_ind = interp1(x(~isnan(pred1)),pred1(~isnan(pred1)),x_interp);
                ca1_ind = x_interp([find(ca1_ind>0.4,1,'first') find(ca1_ind<0.6,1,'last')]);

                subset = ca1gam(:,1)>=ca1_ind(1) & ca1gam(:,1)<=ca1_ind(2);
                rip_subset = times(subset);

                ca1_rip = ca1gam(rip_subset==1,2); ca1_norip = ca1gam(rip_subset==0,2);

                if length(ca1_rip)>20 && length(ca1_norip)>20
                    %Compute the modulation of CA1 and CA3 neurons
                    cellindex = evaluatefilter(cellinfo{day}{epoch},'$meanrate < 7');

                    %Define 200ms time windows
                    riptimes = [ca1_rip ca1_rip+0.2];
                    noriptimes = [ca1_norip ca1_norip+0.2];

                    rip = []; norip = [];
                    for cellcount = 1:size(cellindex,1)
                        ind = [day epoch cellindex(cellcount,:)];

                        if ~isempty(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data)
                            spiketimes = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
                            if isfield(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)},'globalgammaphase')
                                rgam = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.globalgammaphase;
                            else
                                rgam = nan(size(spiketimes));
                            end
                        else
                            spiketimes = [];
                            rgam = [];
                        end
                        %Find valid spikes
                        rip_spikebins = periodAssign(spiketimes, riptimes);
                        norip_spikebins = periodAssign(spiketimes, noriptimes);

                        nrspiketimes = spiketimes; nrgam = rgam;
                        if ~isempty(spiketimes)
                            spiketimes = spiketimes(find(rip_spikebins));
                            rgam = rgam(find(rip_spikebins));
                            nrspiketimes = nrspiketimes(find(norip_spikebins));
                            nrgam = nrgam(find(norip_spikebins));

                            tmpripdata = rgam;
                            tmpnoripdata = nrgam;
                        else
                            tmpripdata = NaN;
                            tmpnoripdata = NaN;
                        end

                        %add cell location: 1 if CA1, 3 if CA3
                        if ~isfield(cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)},'area')
                            tmpripdata(:,2) = 0;
                            tmpnoripdata(:,2) = 0;
                        elseif isequal('CA1',cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)}.area)
                            tmpripdata(:,2) = 1;
                            tmpnoripdata(:,2) = 1;
                        elseif isequal('CA3',cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)}.area)
                            tmpripdata(:,2) = 3;
                            tmpnoripdata(:,2) = 3;
                        else
                            tmpripdata(:,2) = 0;
                            tmpnoripdata(:,2) = 0;
                        end

                        rip = [rip; tmpripdata];
                        norip = [norip; tmpnoripdata];
                    end                    
                    out{day}{epoch}.ripple_ca1_mod = rip(rip(:,2) == 1,1);
                    out{day}{epoch}.ripple_ca3_mod = rip(rip(:,2) == 3,1);
                    out{day}{epoch}.no_ripple_ca1_mod = norip(norip(:,2) == 1,1);
                    out{day}{epoch}.no_ripple_ca3_mod = norip(norip(:,2) == 3,1);

                end
            end
        end
	end
end
end

            