function stats = compute_signal_to_noise(directoryname,fileprefix,days,varargin)

%
%Goes through each day and computes the coherence and phase for equivalent
%powers of ca1 and ca3 gamma 0.5 in the presence and absense of a ripple

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

stats = cell(max(days),1);

%Normalize by the mean and std of the first session
mean_ca1 = []; std_ca1 = [];
days = days(:)';
for day = days
    %load the ripples, tetinfo, and cellinfo files
    tetinfo = loaddatastruct(directoryname,fileprefix,'tetinfo');
    cellinfo = loaddatastruct(directoryname,fileprefix,'cellinfo');
    ripples = loaddatastruct(directoryname,fileprefix,'ripples');
    stats{day} = cell(min(7,length(ripples{day})),1);
    for epoch = 1:min(7,length(ripples{day}))
        
        %initialize output
        stats{day}{epoch}.ripple_ca1_ca3_coherence = [];
        stats{day}{epoch}.ripple_ca1_ca3_phase = [];
        stats{day}{epoch}.no_ripple_ca1_ca3_coherence = [];
        stats{day}{epoch}.no_ripple_ca1_ca3_phase = [];
        
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
       
        %Compute GLM model to determine significance
        [beta dev tmpstats] = glmfit(ca1gam,times,'binomial');

        %Make histogram of data to show where gamma power predicts ripples
        subs = lookup(ca1gam(:,1),x);
        pred1 = accumarray(subs,times,[length(x) 1],@(x) mean(x), NaN);
        
        stats{day}{epoch}.ca1_beta = beta;
        stats{day}{epoch}.ca1_pvalue = tmpstats.p;
        stats{day}{epoch}.ca1_prediction = pred1;
        
        %Go through and pick out 200 ms bins that correspond to 40-60%
        %predictability. Divide into rip+ and rip- times.
        if max(pred1 > 0.6)
            x_interp = -10:0.1:20;
            ca1_ind = interp1(x(~isnan(pred1)),pred1(~isnan(pred1)),x_interp);
            ca1_ind = x_interp([find(ca1_ind>0.4,1,'first') find(ca1_ind<0.6,1,'last')]);
        
            subset = ca1gam(:,1)>=ca1_ind(1) & ca1gam(:,1)<=ca1_ind(2);
            time_subset = time(1:bin:end);
            time_subset = time_subset(subset);
            rip_subset = times(subset);
            
            %Align to peak gamma power
            ca1_rip = ca1gam(rip_subset==1,2); ca1_norip = ca1gam(rip_subset==0,2);

            if sum(rip_subset==1)>0 && sum(rip_subset==0)>0
                [h p] = ttest2(ca1gam(rip_subset==1,1),ca1gam(rip_subset==0,1));
                stats{day}{epoch}.gam_power_rip_norip_pvalue = p;
            else
                stats{day}{epoch}.gam_power_rip_norip_pvalue = NaN;
            end
            if length(ca1_rip)>20 && length(ca1_norip)>20
                %Go through and compute the average ca1-ca1, ca1-ca3, and ca3-ca3
                %coherence and phase for rip and no rip
                params = {};
                params.Fs = 1500;
                params.fpass = [20 50];
                params.trialave =0;
                win = [0 0.1];
                cwin = [0.1 0.01];

                rip_c13 = []; rip_p13 = []; norip_c13 = []; norip_p13 = [];
                min_endtime = time(end);
                for i = ca1tetrodes'
                    loadfile = sprintf('%s/EEGnonreference/%seeg%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,i);
                    load(loadfile)
                    e1eeg = eeg{day}{epoch}{i};
                    e1times = geteegtimes(e1eeg);

                    for j = ca3tetrodes'
                        loadfile = sprintf('%s/EEGnonreference/%seeg%02d-%d-%02d.mat', directoryname, fileprefix,day,epoch,j);
                        load(loadfile)
                        e2 = eeg{day}{epoch}{j};
                        e2times = geteegtimes(e2); clear eeg
                        if length(e1times)>length(e2times)
                            temp = lookup(e2times,e1times);
                            e1times = e1times(temp);
                            e1 = e1eeg.data(temp);
                            e2 = e2.data;
                        elseif length(e2times)>length(e1times)
                            temp = lookup(e1times,e2times);
                            e1 = e1eeg.data;
                            e2 = e2.data(temp);
                        elseif all(isequal(e1times,e2times))
                            e1 = e1eeg.data;
                            e2 = e2.data;
                        else
                            invalid_start = lookup(e1times(1),e2times);
                            if invalid_start > 1
                                e2times = e2times(invalid_start:end);
                                e2.data = e2.data(invalid_start:end);
                            end
                            invalid_start = lookup(e2times(1),e1times);
                            if invalid_start > 1
                                e1times = e1times(invalid_start:end);
                                e1eeg.data = e1eeg.data(invalid_start:end);
                            end
                            invalid_end = lookup(e1times(end),e2times);
                            if ~isequal(invalid_end,length(e2times))
                                e2times = e2times(1:invalid_end);
                                e2.data = e2.data(1:invalid_end);
                            end
                            invalid_end = lookup(e2times(end),e1times);
                            if ~isequal(invalid_end,length(e1times))
                                e1times = e1times(1:invalid_end);
                                e1eeg.data = e1eeg.data(1:invalid_end);
                            end
                            e1 = e1eeg.data;
                            e2 = e2.data;
                        end

                        starttime = e1times(1);
                        
                        %Define triggering events
                        triggers_rip = ca1_rip-starttime;
                        triggers_norip = ca1_norip-starttime;

                        %Remove triggering events that are too close to the beginning or end
                        while triggers_rip(1)<win(1) + 0.5
                            triggers_rip(1) = [];
                        end
                        while triggers_norip(1)<win(1) + 0.5
                            triggers_norip(1) = [];
                        end
                        while triggers_rip(end)> min_endtime-win(2)-0.5
                            triggers_rip(end) = [];
                        end
                        while triggers_norip(end)> min_endtime-win(2)-0.5
                            triggers_norip(end) = [];
                        end

                        % Calculate the event triggered coherence
                        data1 = createdatamatc(e1,triggers_rip,params.Fs,[win(1) win(2)]);
                        data2 = createdatamatc(e2,triggers_rip,params.Fs,[win(1) win(2)]);

                        [C,phi] = coherencyc(data1,data2,params);

                        if isempty(rip_c13)
                            rip_c13 = mean(C);
                            count_rip13 = 1;
                            rip_p13 = mean(phi);
                        else
                            rip_c13 = rip_c13+mean(C);
                            count_rip13 = count_rip13+1;
                            rip_p13 = cat(1,rip_p13,mean(phi));
                        end

                        % Calculate the event triggered coherence
                        data1 = createdatamatc(e1,triggers_norip,params.Fs,[win(1) win(2)]);
                        data2 = createdatamatc(e2,triggers_norip,params.Fs,[win(1) win(2)]);

                        [C,phi] = coherencyc(data1,data2,params);

                        if isempty(norip_c13)
                            norip_c13 = mean(C);
                            count_norip13 = 1;
                            norip_p13 = mean(phi);
                        else
                            norip_c13 = norip_c13+mean(C);
                            count_norip13 = count_norip13+1;
                            norip_p13 = cat(1,norip_p13,mean(phi));
                        end
                    end
                end
                rip_c13 = rip_c13./count_rip13; norip_c13 = norip_c13./count_norip13;
                rip_p13 = mean(rip_p13);    norip_p13 = mean(rip_p13);
                stats{day}{epoch}.ripple_ca1_ca3_coherence = rip_c13;
                stats{day}{epoch}.ripple_ca1_ca3_phase = rip_p13;
                stats{day}{epoch}.no_ripple_ca1_ca3_coherence = norip_c13;
                stats{day}{epoch}.no_ripple_ca1_ca3_phase = norip_p13;
            end            
        end
        end
    end
end

            