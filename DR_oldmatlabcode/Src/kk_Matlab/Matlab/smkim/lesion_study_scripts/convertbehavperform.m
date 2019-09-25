
dayepoch = [ 1 1; ...
    1 2; ...
    2 1; ...
    2 2; ...
    3 1; ...
    3 2; ...
    4 1; ...
    4 2; ...
    5 1; ...
    5 2; ...
    6 1; ...
    6 2; ...
    7 1; ...
    7 2; ...
    8 1; ...
    8 2; ...
    9 1; ...
    9 2; ...
    10 1; ...
    10 2 ];
criterion = 1/2;
conf = 0.95;

load('groups.mat');
i = 0;
for s = 1:length(groups)
    if groups(s).keep
        i = i+1;
        behavperform(i).subject = groups(s).subject;
        behavperform(i).group = groups(s).group;
        behavperform(i).dayepoch = dayepoch;
        behavperform(i).inreward = zeros(0,1);
        behavperform(i).outreward = zeros(0,1);
        behavperform(i).dayintrials = nan([size(dayepoch,1),2]);
        behavperform(i).dayouttrials = nan([size(dayepoch,1),2]);
        behavperform(i).conf = conf;
        behavperform(i).criterion = criterion;
        load([behavperform(i).subject '_Wtrack_journeys.mat']);
        for j = 1:size(dayepoch,1)
            traj = journeys{dayepoch(j,1)}{dayepoch(j,2)};
            firsttrial = NaN; lasttrial = NaN;    
            for k = 1:length(traj.correct)
                if ( (traj.startzone(k)==1) | (traj.startzone(k)==3) ) ...
                    & ~isnan(traj.endzone(k)) & ~isnan(traj.correct(k))
                    behavperform(i).inreward(end+1,1) = traj.correct(k);
                    if isnan(behavperform(i).dayintrials(j,1))
                        behavperform(i).dayintrials(j,1) = numel(behavperform(i).inreward);
                    end
                    behavperform(i).dayintrials(j,2) = numel(behavperform(i).inreward);
                elseif (traj.startzone(k)==2) & ~isnan(traj.endzone(k)) ...
                    & ~isnan(traj.correct(k))
                    behavperform(i).outreward(end+1,1) = traj.correct(k);
                    if isnan(behavperform(i).dayouttrials(j,1))
                        behavperform(i).dayouttrials(j,1) = numel(behavperform(i).outreward);
                    end
                    behavperform(i).dayouttrials(j,2) = numel(behavperform(i).outreward);
                end
            end
        end
        clear('journeys');

        % now estimate probability correct
        behavperform(i).inprobcorrect = getestprobcorrect(behavperform(i).inreward, criterion, 1, conf)
        behavperform(i).outprobcorrect = getestprobcorrect(behavperform(i).outreward, criterion, 1, conf)

  
        % identify learning trial
        lt_in = find(behavperform(s).inprobcorrect(:,2) <= ...
            behavperform(s).criterion,1,'last');
        if (lt_in == size(behavperform(s).inprobcorrect,1));
          behavperform(s).lt_in = -lt_in; 
          behavperform(s).ld_in = NaN;
        else
          behavperform(s).lt_in = lt_in;
          d1 = ceil(find(behavperform(s).dayintrials(:,2) >= ...
              behavperform(s).lt_in,1,'first')/2);
          d2 = ceil(find(behavperform(s).dayintrials(:,1) <= ...
              behavperform(s).lt_in,1,'last')/2);
          if (d1 == d2)
            behavperform(s).ld_in = d1;
          else
            error('something went wrong');
          end
        end

        lt_out = find(behavperform(s).outprobcorrect(:,2) <= ...
            behavperform(s).criterion,1,'last');
        if (lt_out == size(behavperform(s).outprobcorrect,1));
          behavperform(s).lt_out = -lt_out;
          behavperform(s).ld_out = NaN;
        else
          behavperform(s).lt_out = lt_out;
          d1 = ceil(find(behavperform(s).dayouttrials(:,2) >= ...
              behavperform(s).lt_out,1,'first')/2);
          d2 = ceil(find(behavperform(s).dayouttrials(:,1) <= ...
              behavperform(s).lt_out,1,'last')/2);
          if (d1 == d2)
            behavperform(s).ld_out = d1;
          else
            error('something went wrong');
          end
        end

    end
end

