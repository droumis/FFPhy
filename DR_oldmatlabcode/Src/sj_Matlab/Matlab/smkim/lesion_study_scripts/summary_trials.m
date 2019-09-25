
%

subjects= { 'M06' 'M24' 'M25' 'M26' '' ...
    'M02' 'M03' 'M12' 'M13' 'M14' 'M16' 'M17' 'M19' 'M20' 'M22' };

% list subjects along rows, grouping controls (top) from lesion (bottom)

% concatenate total running time for each animal, delimiting sessions with
% little green vertical marks

% correct are o, error are x



% iterate through all journeys in the session and check whether endpoint is
% nonzero (i.e. whether the rat arrived at a food well at the end

tmax = 1;

for i = 1:length(subjects)

    % load journey data for the current subject
    currentsubject = subjects{i};
    if isempty(currentsubject)
        continue;
    end
    currentjourneys = [ currentsubject '_Wtrack_journeys.mat'];
    load(currentjourneys);

    trials(i).subject = currentsubject;
    trials(i).correct = []; % times when the rat arrived at food wells and got reward
    trials(i).error = []; % times when the rat arrived at food wells and did not get reward

    % figure out the durations of the sessions (in minutes)
    session_durations(1) = 0;
    for j = 1:10
        for k = 1:2
            session_durations(2*j+k-1) = 1/60 * ... 
                ( journeys{j}{k}(end).endtime - journeys{j}{k}(1).starttime );
        end
    end
    % and use these to generate offsets for the vector of times
    trials(i).toffsets = cumsum(session_durations);
    if ceil(trials(i).toffsets(end)) > tmax
        tmax = ceil(trials(i).toffsets(end));
    end

    % now accumulate vectors of trial times
    for j = 1:10
        for k = 1:2
            last_side_arm = 0;
            for l = 1:numel(journeys{j}{k})
                % score each journey
                score = NaN; % placeholder initial value
                if journeys{j}{k}(l).startpoint & journeys{j}{k}(l).endpoint
                    if journeys{j}{k}(l).endpoint == journeys{j}{k}(l).startpoint
                        score = 0;
                    else
                        if (journeys{j}{k}(l).startpoint == 1) | (journeys{j}{k}(l).startpoint == 3)
                            % test against the return-to-center rule
                            if journeys{j}{k}(l).endpoint == 2
                                score = 1;
                                last_side_arm = journeys{j}{k}(l).startpoint;
                            else % the rat went from side arm to side arm
                                score = 0;
                                last_side_arm = journeys{j}{k}(l).endpoint;                    
                            end
                        elseif journeys{j}{k}(l).startpoint == 2
                            if journeys{j}{k}(l).endpoint == last_side_arm
                                score = 0;
                            else
                                score = 1;
                            end
                            last_size_arm = journeys{j}{k}(l).endpoint;
                        end
                    end
                    if score
                        trials(i).correct(end+1) = 1/60 * ...
                            (journeys{j}{k}(l).endtime - journeys{j}{k}(1).starttime) ...
                            + trials(i).toffsets(2*j+k-2);
                    else
                        trials(i).error(end+1) = 1/60 * ...
                            (journeys{j}{k}(l).endtime - journeys{j}{k}(1).starttime) ...
                            + trials(i).toffsets(2*j+k-2);
                    end
                end

            end
        end
    end
    % count number of food well visits per minute
    trials(i).errorcounts = histc(trials(i).error,0:1:tmax);
    trials(i).correctcounts = histc(trials(i).correct,0:1:tmax);
    trials(i).totalcounts = trials(i).errorcounts + trials(i).correctcounts;

    clear('journeys');
end
%

figure(1);
set(gcf,'Color','k');
set(gca,'YTick',1:length(subjects),'YTickLabel', subjects, ... 
    'Color','k','XColor','w','YColor','w','LineWidth',2,'FontSize',16, ...
    'XLim',[0 tmax],'YLim',[0.5 length(subjects)+0.5]);
set(get(gca,'XLabel'),'Color','w','FontSize',16, ...
    'String','cumulative experience (minutes)');

for i = 1:length(subjects)

    if isempty(subjects{i})
        continue;
    end
    % draw horizontal rule lines
    line([0 trials(i).toffsets(end)],[i i],'Color','w','LineStyle','-','Marker','none','LineWidth',2);

    % draw color densities to represent task performance
    for l = 1:numel(0:1:tmax)
        if trials(i).totalcounts(l)
            currentcolor = 1/trials(i).totalcounts(l) * ...
                [trials(i).errorcounts(l) trials(i).correctcounts(l) trials(i).errorcounts(l) ];
            line([l l],[i-0.03*trials(i).errorcounts(l) i+0.03*trials(i).correctcounts(l)],'Color',currentcolor,'LineWidth',2);
        end
    end


end
% label groups
text('Position',[-25 2.5 0],'String','sham-operated control',...
    'FontSize',16,'Color','w','Rotation',90,'HorizontalAlignment','center');
text('Position',[-25 10.5 0],'String','lesion',...
    'FontSize',16,'Color','w','Rotation',90,'HorizontalAlignment','center');

%

% determine the smallest number of cumulative trials completed by any individual rat
smallest_total = Inf;
for i = 1:length(subjects)
    if isempty(subjects{i})
        continue
    end
    if (numel(trials(i).error) + numel(trials(i).correct)) < smallest_total
        smallest_total = numel(trials(i).error) + numel(trials(i).correct);
    end
end

smallest_total = 600;

% allocate an imdata array
imdata = zeros(3*length(subjects),smallest_total,3);

figure(2);

% now plot color-coded cells for correct versus incorrect trials
for i = 1:length(subjects)
    if isempty(subjects{i})
        continue;
    end
    trialsort = zeros(numel(trials(i).correct)+numel(trials(i).error),2);
    trialsort(:,1) = [trials(i).correct'; trials(i).error'];
    trialsort(1:numel(trials(i).correct),2) = 1;
    trialsort = sortrows(trialsort,1);
    trialsort = trialsort(1:smallest_total,:);
    imdata(3*i-1,find(trialsort(:,2)),2) = 1;
    imdata(3*i-1,find(~trialsort(:,2)),1) = 1;
    imdata(3*i-1,find(~trialsort(:,2)),3) = 1;
end

image([1 smallest_total],[1-1/3 length(subjects)+1/3],imdata);
set(gcf,'Color','k');
set(gca,'YTick',1:length(subjects),'YTickLabel', subjects, ... 
    'Color','k','XColor','w','YColor','w','LineWidth',2,'FontSize',16, ...
    'XLim',[0 smallest_total+1],'YLim',[0.5 length(subjects)+0.5], ...
    'Box','off','YDir','normal');
set(get(gca,'XLabel'),'Color','w','FontSize',16, ...
    'String','cumulative trials');

% label groups
text('Position',[-50 2.5 0],'String','sham-operated control',...
    'FontSize',16,'Color','w','Rotation',90,'HorizontalAlignment','center');
text('Position',[-50 10.5 0],'String','lesion',...
    'FontSize',16,'Color','w','Rotation',90,'HorizontalAlignment','center');

%}
